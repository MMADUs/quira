#include <iostream>
#include <quira/quira.hpp>

int main() {
  quira::QuantumCircuit circuit(3, 3);

  // Qubit layout:
  //   q0 = message qubit
  //   q1 = Alice's entangled qubit
  //   q2 = Bob's entangled qubit
  //
  // Prepare |+> on q0. After teleportation, q2 should also be |+>.
  circuit.add(quira::Hadamard(0));

  // Create Bell pair between q1 and q2.
  circuit.add(quira::Hadamard(1));
  circuit.add(quira::ControlledNot(1, 2));

  // Bell-basis measurement between message qubit q0 and Alice qubit q1.
  circuit.add(quira::ControlledNot(0, 1));
  circuit.add(quira::Hadamard(0));
  circuit.measure(0, 0);
  circuit.measure(1, 1);

  // Classical feed-forward corrections on Bob's qubit.
  circuit.conditional_add(1, true, quira::PauliX(2));
  circuit.conditional_add(0, true, quira::PauliZ(2));

  // Verify the teleported |+> state by measuring Bob's qubit in the X basis.
  // Measuring |+> in the X basis should produce 0.
  circuit.add(quira::Hadamard(2));
  circuit.measure(2, 2);

  quira::StateVectorSimulator simulator(42);

  const quira::SimulationOutput output = simulator.run(circuit, 1024);

  std::cout << "Teleportation measurement counts\n";
  std::cout << "Expected: Bob's verification bit should be 0.\n\n";

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
