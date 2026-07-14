#include <iostream>
#include <quira/quira.hpp>

int main() {
  quira::QuantumCircuit circuit(2, 2);

  circuit.add(quira::Hadamard(0));
  circuit.add(quira::ControlledNot(0, 1));
  circuit.measure_all();

  quira::StateVectorSimulator simulator(42);

  const quira::SimulationOutput output = simulator.run(circuit, 1024);

  std::cout << "Bell state measurement counts\n";

  for (const auto& [bit_string, count] : output.get_counts()) {
    std::cout << bit_string << ": " << count << '\n';
  }

  const auto most_frequent = output.get_most_frequent();

  if (most_frequent) {
    std::cout << "Most frequent: |" << most_frequent->first
              << ">: " << most_frequent->second << '\n';
  }

  return 0;
}
