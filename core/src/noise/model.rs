use rand::seq::IndexedRandom;

use crate::{bit::QuantumBit, io::QuantumState, ops::QuantumGate};

pub struct NoiseChannel<T: NoiseChannelOperation> {
    channel: T,
    gates: Option<Vec<String>>,
    selective: Option<Vec<QuantumBit>>,
    idle_delay: Option<usize>,
    is_configured: bool,
}

impl<T: NoiseChannelOperation> NoiseChannel<T> {
    pub fn new(channel: T) -> Self {
        Self {
            channel,
            gates: None,
            selective: None,
            idle_delay: None,
            is_configured: false,
        }
    }

    pub fn global_error(&mut self) {
        self.selective = None;
        self.is_configured = true;
    }

    pub fn selective_error<I>(&mut self, qubits: I)
    where
        I: IntoIterator<Item = QuantumBit>,
    {
        self.selective = Some(qubits.into_iter().collect());
        self.is_configured = true;
    }

    pub fn specified_gates<I>(&mut self, gates: I)
    where
        I: IntoIterator<Item = String>,
    {
        let added_gates: Vec<String> = gates.into_iter().collect();
        if let Some(existing_gates) = self.gates.as_mut() {
            existing_gates.extend(added_gates);
        } else {
            self.gates = Some(added_gates);
        }
    }

    pub fn set_delay(&mut self, delay: usize) {
        self.idle_delay = Some(delay);
    }
}

// TODO: 
// 1. preset the noise configuration on each of the channel itself (Delete noise channnel struct)
// 2. build it onwards for the noise channel itself 
// 3. since density is the only support i knew,  make a trait for each qstate for compatibility
// functions.
pub struct NoiseModel<T: NoiseChannelOperation> {
    channels: Vec<NoiseChannel<T>>,
    num_qubits: usize,
}

impl<T: NoiseChannelOperation> NoiseModel<T> {
    pub fn new() -> Self {
        Self {
            channels: Vec::new(),
            num_qubits: 0,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.channels.is_empty()
    }

    pub fn extend(&mut self, model: NoiseModel<T>) {
        self.channels.extend(model.channels);
    }

    pub fn add_error<QI, GI>(&mut self, error: T, qubits: QI, gates: GI, delay: Option<usize>)
    where
        QI: IntoIterator<Item = QuantumBit>,
        GI: IntoIterator<Item = String>,
    {
        let mut channel = NoiseChannel::new(error);

        let qubit_vec: Vec<QuantumBit> = qubits.into_iter().collect();
        let gate_vec: Vec<String> = gates.into_iter().collect();

        if qubit_vec.is_empty() {
            channel.global_error();
        } else {
            channel.selective_error(qubit_vec);
        }

        if !gate_vec.is_empty() {
            channel.specified_gates(gate_vec);
        }

        if let Some(d) = delay {
            channel.set_delay(d);
        }

        self.channels.push(channel);
    }

    pub(crate) fn build_noise(&mut self, n: usize) {
        self.num_qubits = n;
        // build all the noise channel in the noise model
    }

    pub(crate) fn simulate_noise<Q: QuantumState>(
        &self,
        state: &mut Q,
        ops: &Box<dyn QuantumGate>,
    ) {
        // 1. check for construct targets, make it only support 1 and 2 qubit gate.
        // 2. check for state, make it only support for density.
        // 3. validate if the noise channel is applied for specified gates.
        // 4. perform the noise model.
    }
}

pub trait NoiseChannelOperation {
    fn is_idle_compatible(&self) -> bool;

    fn global(&self, num_qubits: usize) -> Self;

    fn selective(&self, total_qubits: usize, specified: &[QuantumBit]) -> Self;
}

#[derive(Clone, Debug)]
pub struct NoNoise;

impl NoiseChannelOperation for NoNoise {
    fn is_idle_compatible(&self) -> bool {
        false
    }

    fn global(&self, _num_qubits: usize) -> Self {
        Self
    }

    fn selective(&self, _total_qubits: usize, _specified: &[QuantumBit]) -> Self {
        Self
    }
}
