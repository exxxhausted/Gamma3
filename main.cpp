#include <iostream>
#include <random>
#include <numbers>

#include <cxxopts.hpp>

#include "Gamma3.hpp"

using namespace gamma3;

#define N_A 6.02214076e23
#define RAY_COUNT 100000
#define DISTANCE 100
#define ENERGY_MEV 0.662

physics::Photon spawn_photon_at_point(const geometry::Point& p, double E_MeV);
geometry::Surface create_cube_detector_mm(double cube_size_mm = 100.0, const geometry::Point& center = geometry::Point(0, 0, 0));
double compute_sigma_Atom(double E_MeV, const materials::Atom& atom, etc::Function2<double>& f);
double compute_sigma_Molecule(double E_MeV, const materials::Molecule& molecule, etc::Function2<double>& f);
double compute_sigma_Material(double E_MeV, const materials::Material& material, etc::Function2<double>& f);

int main(int argc, char** argv) {

    //CLI INTERFACE

    cxxopts::Options options("Gamma3", "Gamma simulation");

    options.add_options()
        ("r,repo", "Path to Database.json", cxxopts::value<std::string>()->default_value("../../Database.json"))
        ("o,output", "Output histogram file", cxxopts::value<std::string>()->default_value("histogram.txt"))
        ("h,help", "Print help");

    auto result = options.parse(argc, argv);

    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        return 0;
    }

    std::string repo_path = result["repo"].as<std::string>();
    std::string output_path = result["output"].as<std::string>();


    //SIMULATION
    try {
        auto start = std::chrono::high_resolution_clock::now();

        auto repo = io::JsonRepository(repo_path);

        auto NaI = repo.loadMaterial("NaI");
        auto sigma_photo_formula = repo.loadFormula2("SigmaPhoto");
        auto sigma_compton_formula = repo.loadFormula2("SigmaCompton");
        auto sigma_pair_formula = repo.loadFormula2("SigmaPair");

        auto cube = create_cube_detector_mm();

        auto gen = std::mt19937(std::random_device{}());
        auto dist = std::uniform_real_distribution<double>(0, 1);
        int cnt = 0;
        int absorbed_photons = 0;
        std::vector<double> histohram;

        for (int i = 0; i < RAY_COUNT; ++i) {
            auto ph = spawn_photon_at_point(geometry::Point(50 + DISTANCE, 0, 0), ENERGY_MEV);
            auto intersection = cube.intersect(ph.ray());
            if (intersection != std::nullopt) {
                cnt++;
                ph.move(intersection->t + 0.0000001);

                bool photon_absorbed = false;

                double delta_E = 0.0;

                while (cube.contains(ph.ray().source()) && !photon_absorbed) {
                    auto gamma = dist(gen);

                    double sigma_photo = compute_sigma_Material(ph.energy(), NaI, sigma_photo_formula);
                    double sigma_compton = compute_sigma_Material(ph.energy(), NaI, sigma_compton_formula);
                    double sigma_pair = compute_sigma_Material(ph.energy(), NaI, sigma_pair_formula);
                    double sigma = sigma_photo + sigma_compton + sigma_pair;

                    double lambda = -(1 / sigma) * log(gamma);
                    ph.move(lambda);

                    if (!cube.contains(ph.ray().source())) break;

                    double P_photo = sigma_photo / sigma;
                    double P_compton = sigma_compton / sigma;
                    double P_pair = sigma_pair / sigma;
                    double random_value = dist(gen);

                    if (random_value < P_photo) {
                        // Photoabsorption
                        delta_E += ph.energy();
                        absorbed_photons++;
                        photon_absorbed = true;
                    }
                    else if (random_value < P_photo + P_compton) {
                        // Compton scattering
                        double alpha = ph.energy() / 0.511;
                        auto old_dir = glm::normalize(ph.ray().direction());
                        auto new_dir = glm::normalize(geometry::Vector(dist(gen), dist(gen), dist(gen)));
                        double cos_theta = old_dir.x * new_dir.x + old_dir.y * new_dir.y + old_dir.z * new_dir.z;
                        double new_energy = ph.energy() / (1.0 + alpha * (1.0 - cos_theta));
                        double energy_loss = ph.energy() - new_energy;

                        ph.setDirection(new_dir);
                        ph.setEnergy(new_energy);
                        delta_E += energy_loss;
                    }
                    else if (random_value < P_photo + P_compton + P_pair) {
                        // Pair production
                        delta_E += ph.energy();
                        absorbed_photons++;
                        photon_absorbed = true;
                    }
                }
                if (delta_E != 0.0) histohram.push_back(delta_E);
            }
        }

        std::cout << "P_absorption: " << static_cast<double>(absorbed_photons) / cnt << std::endl;

        std::cout << "Intersection counter: " << cnt << std::endl;

        {
            std::ofstream fout(output_path);
            if (!fout) {
                std::cerr << "Cannot open output file: " << output_path << std::endl;
            } else {
                for (double v : histohram) fout << v << "\n";
                std::cout << "Histogram saved to: " << output_path << std::endl;
            }
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Execution time: " << duration.count() << " ms" << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }
    return 0;
}

// HELPER FUNCTIONS

physics::Photon spawn_photon_at_point(const geometry::Point& p, double E_MeV) {
    std::random_device rd;
    static std::mt19937 rng(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    double theta = 2.0 * std::numbers::pi * dist(rng);
    double phi = std::acos(1.0 - 2.0 * dist(rng));

    geometry::Vector direction(
        std::sin(phi) * std::cos(theta),
        std::sin(phi) * std::sin(theta),
        std::cos(phi)
        );

    return physics::Photon(geometry::Ray(p, direction), E_MeV);
}

geometry::Surface create_cube_detector_mm(double cube_size_mm, const geometry::Point& center) {
    double half_size = cube_size_mm / 2.0;

    std::vector<geometry::Point> vertices = {
        {center.x - half_size, center.y - half_size, center.z - half_size},
        {center.x + half_size, center.y - half_size, center.z - half_size},
        {center.x + half_size, center.y + half_size, center.z - half_size},
        {center.x - half_size, center.y + half_size, center.z - half_size},

        {center.x - half_size, center.y - half_size, center.z + half_size},
        {center.x + half_size, center.y - half_size, center.z + half_size},
        {center.x + half_size, center.y + half_size, center.z + half_size},
        {center.x - half_size, center.y + half_size, center.z + half_size}
    };

    std::vector<glm::u32vec3> faces = {
        {0, 2, 1}, {0, 3, 2},
        {4, 5, 6}, {4, 6, 7},
        {0, 1, 5}, {0, 5, 4},
        {1, 2, 6}, {1, 6, 5},
        {2, 3, 7}, {2, 7, 6},
        {3, 0, 4}, {3, 4, 7}
    };

    return geometry::Surface(vertices, faces);
}

double compute_sigma_Atom(double E_MeV,
                          const materials::Atom& atom,
                          etc::Function2<double>& f) {
    return f(atom.Z(), E_MeV);
}

double compute_sigma_Molecule(double E_MeV,
                              const materials::Molecule& molecule,
                              etc::Function2<double>& f) {
    auto& composition = molecule.components();
    double sum = 0.0;
    for(const auto& [atom, count] : composition) sum += count * compute_sigma_Atom(E_MeV, atom, f);
    return sum;
}

double compute_sigma_Material(double E_MeV,
                              const materials::Material& material,
                              etc::Function2<double>& f) {
    auto& composition = material.components();
    double sum = 0.0;
    for(const auto& [molecule, mass_fraction] : composition) {
        double M = molecule.mass();
        sum += mass_fraction * (N_A / M) * compute_sigma_Molecule(E_MeV, molecule, f);
    }
    return sum * material.density();
}

