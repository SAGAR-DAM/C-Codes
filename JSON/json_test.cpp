#include <iostream>
#include <fstream>
#include "json.hpp"

using json = nlohmann::json;

int main() {
    std::ifstream file("data.json");
    if (!file.is_open()) {
        std::cerr << "Cannot open JSON file.\n";
        return 1;
    }

    json j;
    file >> j;

    std::cout << "Name: " << j["name"] << "\n";
    std::cout << "Age: " << j["age"] << "\n";

    std::cout << "Skills: ";
    for (const auto& skill : j["skills"]) {
        std::cout << skill << " ";
    }
    std::cout << "\n";

    std::cout << "Institution: " << j["details"]["institution"] << "\n";
    std::cout << "Project: " << j["details"]["project"] << "\n";

    return 0;
}
