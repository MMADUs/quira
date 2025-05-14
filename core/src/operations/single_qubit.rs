use super::{QuantumGate, SingleQubit};
use crate::math::{Complex, Matrix, PI};
use ndarray::array;

pub enum SingleQubitOperation {
    SingleQubitGate,
    RotateX,
    RotateY,
    RotateZ,
    RotateXY,
    PauliX,
    PauliY,
    PauliZ,
    SqrtPauliX,
    InvSqrtPauliX,
    SqrtPauliY,
    InvSqrtPauliY,
    Hadamard,
    SGate,
    InvSGate,
    SXGate,
    InvSXGate,
    TGate,
    InvTGate,
    UGate,
    PhaseShift0,
    PhaseShift1,
    RotateAroundSphericalAxis,
    GPi,
    GPi2,
    Identity,
}

#[derive(Debug, Clone)]
pub struct SingleQubitGate {
    qubit: usize,
    alpha_re: f64,
    alpha_im: f64,
    beta_re: f64,
    beta_im: f64,
    global_phase: f64,
}

impl SingleQubitGate {
    pub fn new(
        qubit: usize,
        alpha_re: f64,
        alpha_im: f64,
        beta_re: f64,
        beta_im: f64,
        global_phase: f64,
    ) -> Self {
        let norm: f64 = alpha_re.powi(2) + alpha_im.powi(2) + beta_re.powi(2) + beta_im.powi(2);
        if (norm - 1.0).abs() > 1e-6 {
            panic!("unitary normalization error");
        }
        Self {
            qubit,
            alpha_re,
            alpha_im,
            beta_re,
            beta_im,
            global_phase,
        }
    }
}

impl QuantumGate for SingleQubitGate {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let pref = Complex::from_polar(1.0, self.global_phase);
        array![
            [
                pref * Complex::new(self.alpha_re, self.alpha_im),
                pref * Complex::new(-1.0 * self.beta_re, self.beta_im)
            ],
            [
                pref * Complex::new(self.beta_re, self.beta_im),
                pref * Complex::new(self.alpha_re, -1.0 * self.alpha_im)
            ]
        ]
    }

    fn name(&self) -> String {
        String::from("Single-Qubit")
    }
}

impl SingleQubit for SingleQubitGate {
    fn target_qubit(&self) -> usize {
        self.qubit
    }
}

#[derive(Debug, Clone)]
pub struct RotateX {
    qubit: usize,
    theta: f64,
}

impl RotateX {
    pub fn new(qubit: usize, theta: f64) -> Self {
        Self { qubit, theta }
    }
}

impl QuantumGate for RotateX {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let c: f64 = (self.theta / 2.0).cos();
        let s: f64 = (self.theta / 2.0).sin();
        array![
            [Complex::new(c, 0.0), Complex::new(0.0, -1.0 * s)],
            [Complex::new(0.0, -1.0 * s), Complex::new(c, 0.0)]
        ]
    }

    fn name(&self) -> String {
        String::from("RX")
    }
}

impl SingleQubit for RotateX {
    fn target_qubit(&self) -> usize {
        self.qubit
    }
}

#[derive(Debug, Clone)]
pub struct RotateY {
    qubit: usize,
    theta: f64,
}

impl RotateY {
    pub fn new(qubit: usize, theta: f64) -> Self {
        Self { qubit, theta }
    }
}

impl QuantumGate for RotateY {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let c: f64 = (self.theta / 2.0).cos();
        let s: f64 = (self.theta / 2.0).sin();
        array![
            [Complex::new(c, 0.0), Complex::new(-1.0 * s, 0.0)],
            [Complex::new(s, 0.0), Complex::new(c, 0.0)]
        ]
    }

    fn name(&self) -> String {
        String::from("RY")
    }
}

impl SingleQubit for RotateY {
    fn target_qubit(&self) -> usize {
        self.qubit
    }
}

#[derive(Debug, Clone)]
pub struct RotateZ {
    qubit: usize,
    theta: f64,
}

impl RotateZ {
    pub fn new(qubit: usize, theta: f64) -> Self {
        RotateZ { qubit, theta }
    }
}

impl QuantumGate for RotateZ {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let c: f64 = (self.theta / 2.0).cos();
        let s: f64 = (self.theta / 2.0).sin();
        array![
            [Complex::new(c, -1.0 * s), Complex::new(0.0, 0.0)],
            [Complex::new(0.0, 0.0), Complex::new(c, s)]
        ]
    }

    fn name(&self) -> String {
        String::from("RZ")
    }
}

impl SingleQubit for RotateZ {
    fn target_qubit(&self) -> usize {
        self.qubit
    }
}

#[derive(Debug, Clone)]
pub struct RotateXY {
    qubit: usize,
    theta: f64,
    phi: f64,
}

impl RotateXY {
    pub fn new(qubit: usize, theta: f64, phi: f64) -> Self {
        Self { qubit, theta, phi }
    }
}

impl QuantumGate for RotateXY {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let c: f64 = (self.theta / 2.0).cos();
        let s: f64 = (self.theta / 2.0).sin();
        let vx: f64 = (self.phi).cos();
        let vy: f64 = (self.phi).sin();
        array![
            [
                Complex::new(c, 0.0),
                Complex::new(-1.0 * s * vy, -1.0 * s * vx)
            ],
            [Complex::new(s * vy, -1.0 * s * vx), Complex::new(c, 0.0)]
        ]
    }

    fn name(&self) -> String {
        String::from("RXY")
    }
}

impl SingleQubit for RotateXY {
    fn target_qubit(&self) -> usize {
        self.qubit
    }
}

#[derive(Debug, Clone)]
pub struct PauliX {
    qubit: usize,
}

impl PauliX {
    pub fn new(qubit: usize) -> Self {
        Self { qubit }
    }
}

impl QuantumGate for PauliX {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        array![
            [Complex::new(0.0, 0.0), Complex::new(1.0, 0.0)],
            [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)]
        ]
    }

    fn name(&self) -> String {
        String::from("X")
    }
}

impl SingleQubit for PauliX {
    fn target_qubit(&self) -> usize {
        self.qubit
    }
}

#[derive(Debug, Clone)]
pub struct PauliY {
    qubit: usize,
}

impl PauliY {
    pub fn new(qubit: usize) -> Self {
        Self { qubit }
    }
}

impl QuantumGate for PauliY {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        array![
            [Complex::new(0.0, 0.0), Complex::new(0.0, -1.0)],
            [Complex::new(0.0, 1.0), Complex::new(0.0, 0.0)]
        ]
    }

    fn name(&self) -> String {
        String::from("Y")
    }
}

impl SingleQubit for PauliY {
    fn target_qubit(&self) -> usize {
        self.qubit
    }
}

#[derive(Debug, Clone)]
pub struct PauliZ {
    qubit: usize,
}

impl PauliZ {
    pub fn new(qubit: usize) -> Self {
        Self { qubit }
    }
}

impl QuantumGate for PauliZ {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        array![
            [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
            [Complex::new(0.0, 0.0), Complex::new(-1.0, 0.0)]
        ]
    }

    fn name(&self) -> String {
        String::from("Z")
    }
}

impl SingleQubit for PauliZ {
    fn target_qubit(&self) -> usize {
        self.qubit
    }
}

#[derive(Debug, Clone)]
pub struct SqrtPauliX {
    qubit: usize,
}

impl SqrtPauliX {
    pub fn new(qubit: usize) -> Self {
        Self { qubit }
    }
}

impl QuantumGate for SqrtPauliX {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let theta: f64 = PI / 2.0;
        let c: f64 = (theta / 2.0).cos();
        let s: f64 = (theta / 2.0).sin();
        array![
            [Complex::new(c, 0.0), Complex::new(0.0, -1.0 * s)],
            [Complex::new(0.0, -1.0 * s), Complex::new(c, 0.0)]
        ]
    }

    fn name(&self) -> String {
        String::from("Sqrt-X")
    }
}

impl SingleQubit for SqrtPauliX {
    fn target_qubit(&self) -> usize {
        self.qubit
    }
}

#[derive(Debug, Clone)]
pub struct InvSqrtPauliX {
    qubit: usize,
}

impl InvSqrtPauliX {
    pub fn new(qubit: usize) -> Self {
        Self { qubit }
    }
}

impl QuantumGate for InvSqrtPauliX {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let theta: f64 = PI / 2.0;
        let c: f64 = (theta / 2.0).cos();
        let s: f64 = (theta / 2.0).sin();
        array![
            [Complex::new(c, 0.0), Complex::new(0.0, 1.0 * s)],
            [Complex::new(0.0, 1.0 * s), Complex::new(c, 0.0)]
        ]
    }

    fn name(&self) -> String {
        String::from("Inv-Sqrt-X")
    }
}

impl SingleQubit for InvSqrtPauliX {
    fn target_qubit(&self) -> usize {
        self.qubit
    }
}

#[derive(Debug, Clone)]
pub struct SqrtPauliY {
    qubit: usize,
}

impl SqrtPauliY {
    pub fn new(qubit: usize) -> Self {
        Self { qubit }
    }
}

impl QuantumGate for SqrtPauliY {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let theta: f64 = PI / 2.0;
        let c: f64 = (theta / 2.0).cos();
        let s: f64 = (theta / 2.0).sin();
        array![
            [Complex::new(c, 0.0), Complex::new(-1.0 * s, 0.0)],
            [Complex::new(s, 0.0), Complex::new(c, 0.0)]
        ]
    }

    fn name(&self) -> String {
        String::from("Sqrt-Y")
    }
}

impl SingleQubit for SqrtPauliY {
    fn target_qubit(&self) -> usize {
        self.qubit
    }
}

#[derive(Debug, Clone)]
pub struct InvSqrtPauliY {
    qubit: usize,
}

impl InvSqrtPauliY {
    pub fn new(qubit: usize) -> Self {
        Self { qubit }
    }
}

impl QuantumGate for InvSqrtPauliY {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let theta: f64 = (-1.0 * PI) / 2.0;
        let c: f64 = (theta / 2.0).cos();
        let s: f64 = (theta / 2.0).sin();
        array![
            [Complex::new(c, 0.0), Complex::new(-1.0 * s, 0.0)],
            [Complex::new(s, 0.0), Complex::new(c, 0.0)]
        ]
    }

    fn name(&self) -> String {
        String::from("Inv-Sqrt-Y")
    }
}

impl SingleQubit for InvSqrtPauliY {
    fn target_qubit(&self) -> usize {
        self.qubit
    }
}

#[derive(Debug, Clone)]
pub struct Hadamard {
    qubit: usize,
}

impl Hadamard {
    pub fn new(qubit: usize) -> Self {
        Self { qubit }
    }
}

impl QuantumGate for Hadamard {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let f = 1.0 / ((2.0_f64).sqrt());
        array![
            [Complex::new(f, 0.0), Complex::new(f, 0.0)],
            [Complex::new(f, 0.0), Complex::new(-1.0 * f, 0.0)]
        ]
    }

    fn name(&self) -> String {
        String::from("H")
    }
}

impl SingleQubit for Hadamard {
    fn target_qubit(&self) -> usize {
        self.qubit
    }
}

#[derive(Debug, Clone)]
pub struct SGate {
    qubit: usize,
}

impl SGate {
    pub fn new(qubit: usize) -> Self {
        Self { qubit }
    }
}

impl QuantumGate for SGate {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        array![
            [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
            [Complex::new(0.0, 0.0), Complex::new(0.0, 1.0)]
        ]
    }

    fn name(&self) -> String {
        String::from("S")
    }
}

impl SingleQubit for SGate {
    fn target_qubit(&self) -> usize {
        self.qubit
    }
}

#[derive(Debug, Clone)]
pub struct InvSGate {
    qubit: usize,
}

impl InvSGate {
    pub fn new(qubit: usize) -> Self {
        Self { qubit }
    }
}

impl QuantumGate for InvSGate {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        array![
            [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
            [Complex::new(0.0, 0.0), Complex::new(0.0, -1.0)]
        ]
    }

    fn name(&self) -> String {
        String::from("Inv-S")
    }
}

impl SingleQubit for InvSGate {
    fn target_qubit(&self) -> usize {
        self.qubit
    }
}

#[derive(Debug, Clone)]
pub struct SXGate {
    qubit: usize,
}

impl SXGate {
    pub fn new(qubit: usize) -> Self {
        Self { qubit }
    }
}

impl QuantumGate for SXGate {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let theta: f64 = PI / 2.0;
        let c: f64 = (theta / 2.0).cos();
        let s: f64 = (theta / 2.0).sin();
        let gp = Complex::new(0.0, PI / 4.0).exp();
        array![
            [gp * Complex::new(c, 0.0), gp * Complex::new(0.0, -1.0 * s)],
            [gp * Complex::new(0.0, -1.0 * s), gp * Complex::new(c, 0.0)]
        ]
    }

    fn name(&self) -> String {
        String::from("SX")
    }
}

impl SingleQubit for SXGate {
    fn target_qubit(&self) -> usize {
        self.qubit
    }
}

pub struct InvSXGate {
    qubit: usize,
}

impl InvSXGate {
    pub fn new(qubit: usize) -> Self {
        Self { qubit }
    }
}

impl QuantumGate for InvSXGate {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let theta: f64 = PI / 2.0;
        let c: f64 = (theta / 2.0).cos();
        let s: f64 = (theta / 2.0).sin();
        let gp = Complex::new(0.0, PI / 4.0).exp();
        array![
            [gp * Complex::new(c, 0.0), gp * Complex::new(0.0, 1.0 * s)],
            [gp * Complex::new(0.0, 1.0 * s), gp * Complex::new(c, 0.0)]
        ]
    }

    fn name(&self) -> String {
        String::from("Inv-SX")
    }
}

impl SingleQubit for InvSXGate {
    fn target_qubit(&self) -> usize {
        self.qubit
    }
}

#[derive(Debug, Clone)]
pub struct TGate {
    qubit: usize,
}

impl TGate {
    pub fn new(qubit: usize) -> Self {
        Self { qubit }
    }
}

impl QuantumGate for TGate {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        array![
            [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
            [
                Complex::new(0.0, 0.0),
                Complex::new((PI / 4.0).cos(), (PI / 4.0).sin())
            ]
        ]
    }

    fn name(&self) -> String {
        String::from("T")
    }
}

impl SingleQubit for TGate {
    fn target_qubit(&self) -> usize {
        self.qubit
    }
}

#[derive(Debug, Clone)]
pub struct InvTGate {
    qubit: usize,
}

impl InvTGate {
    pub fn new(qubit: usize) -> Self {
        Self { qubit }
    }
}

impl QuantumGate for InvTGate {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        array![
            [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
            [
                Complex::new(0.0, 0.0),
                Complex::new((PI / 4.0).cos(), -1.0 * (PI / 4.0).sin())
            ]
        ]
    }

    fn name(&self) -> String {
        String::from("Inv-T")
    }
}

impl SingleQubit for InvTGate {
    fn target_qubit(&self) -> usize {
        self.qubit
    }
}

#[derive(Debug, Clone)]
pub struct UGate {
    qubit: usize,
    theta: f64,
    phi: f64,
    lambda: f64,
}

impl UGate {
    pub fn new(qubit: usize, theta: f64, phi: f64, lambda: f64) -> Self {
        Self {
            qubit,
            theta,
            phi,
            lambda,
        }
    }
}

impl QuantumGate for UGate {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let c: f64 = (self.theta / 2.0).cos();
        let s: f64 = (self.theta / 2.0).sin();

        array![
            [
                Complex::new(c, 0.0),
                Complex::new(-1.0 * s * self.lambda.cos(), -1.0 * s * self.lambda.sin())
            ],
            [
                Complex::new(s * self.phi.cos(), s * self.phi.sin()),
                Complex::new(
                    c * self.phi.cos() * self.lambda.cos() - self.phi.sin() * self.lambda.sin(),
                    c * self.phi.cos() * self.lambda.sin() + self.phi.sin() * self.lambda.cos()
                )
            ]
        ]
    }

    fn name(&self) -> String {
        format!("U({:.4}, {:.4}, {:.4})", self.theta, self.phi, self.lambda)
    }
}

impl SingleQubit for UGate {
    fn target_qubit(&self) -> usize {
        self.qubit
    }
}

#[derive(Debug, Clone)]
pub struct PhaseShiftState1 {
    qubit: usize,
    theta: f64,
}

impl PhaseShiftState1 {
    pub fn new(qubit: usize, theta: f64) -> Self {
        Self { qubit, theta }
    }
}

impl QuantumGate for PhaseShiftState1 {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let theta: f64 = self.theta;
        array![
            [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
            [
                Complex::new(0.0, 0.0),
                Complex::new(theta.cos(), theta.sin())
            ]
        ]
    }

    fn name(&self) -> String {
        String::from("Phase-Shift-1")
    }
}

impl SingleQubit for PhaseShiftState1 {
    fn target_qubit(&self) -> usize {
        self.qubit
    }
}

#[derive(Debug, Clone)]
pub struct PhaseShiftState0 {
    qubit: usize,
    theta: f64,
}

impl PhaseShiftState0 {
    pub fn new(qubit: usize, theta: f64) -> Self {
        Self { qubit, theta }
    }
}

impl QuantumGate for PhaseShiftState0 {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let theta: f64 = self.theta;
        array![
            [
                Complex::new(theta.cos(), theta.sin()),
                Complex::new(0.0, 0.0)
            ],
            [Complex::new(0.0, 0.0), Complex::new(1.0, 0.0)]
        ]
    }

    fn name(&self) -> String {
        String::from("Phase-Shift-0")
    }
}

impl SingleQubit for PhaseShiftState0 {
    fn target_qubit(&self) -> usize {
        self.qubit
    }
}

#[derive(Debug, Clone)]
pub struct RotateAroundSphericalAxis {
    qubit: usize,
    theta: f64,
    spherical_theta: f64,
    spherical_phi: f64,
}

impl RotateAroundSphericalAxis {
    pub fn new(qubit: usize, theta: f64, spherical_theta: f64, spherical_phi: f64) -> Self {
        Self {
            qubit,
            theta,
            spherical_theta,
            spherical_phi,
        }
    }
}

impl QuantumGate for RotateAroundSphericalAxis {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let c: f64 = (self.theta / 2.0).cos();
        let s: f64 = (self.theta / 2.0).sin();
        let vx: f64 = (self.spherical_theta).sin() * (self.spherical_phi).cos();
        let vy: f64 = (self.spherical_theta).sin() * (self.spherical_phi).sin();
        let vz: f64 = (self.spherical_theta).cos();
        array![
            [
                Complex::new(c, -1.0 * s * vz),
                Complex::new(-1.0 * s * vy, -1.0 * s * vx)
            ],
            [Complex::new(s * vy, -1.0 * s * vx), Complex::new(c, s * vz)]
        ]
    }

    fn name(&self) -> String {
        String::from("Rotate-Spherical")
    }
}

impl SingleQubit for RotateAroundSphericalAxis {
    fn target_qubit(&self) -> usize {
        self.qubit
    }
}

#[derive(Debug, Clone)]
pub struct GPi {
    qubit: usize,
    theta: f64,
}

impl GPi {
    pub fn new(qubit: usize, theta: f64) -> Self {
        Self { qubit, theta }
    }
}

impl QuantumGate for GPi {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let c: f64 = (self.theta).cos();
        let s: f64 = (self.theta).sin();
        array![
            [Complex::new(0.0, 0.0), Complex::new(c, -1.0 * s)],
            [Complex::new(c, s), Complex::new(0.0, 0.0)]
        ]
    }

    fn name(&self) -> String {
        String::from("GPi")
    }
}

impl SingleQubit for GPi {
    fn target_qubit(&self) -> usize {
        self.qubit
    }
}

#[derive(Debug, Clone)]
pub struct GPi2 {
    qubit: usize,
    theta: f64,
}

impl GPi2 {
    pub fn new(qubit: usize, theta: f64) -> Self {
        Self { qubit, theta }
    }
}

impl QuantumGate for GPi2 {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let c: f64 = (self.theta).cos();
        let s: f64 = (self.theta).sin();
        array![
            [Complex::new(1.0, 0.0), Complex::new(-1.0 * s, -1.0 * c)],
            [Complex::new(s, -1.0 * c), Complex::new(1.0, 0.0)]
        ] / 2.0_f64.sqrt()
    }

    fn name(&self) -> String {
        String::from("GPi2")
    }
}

impl SingleQubit for GPi2 {
    fn target_qubit(&self) -> usize {
        self.qubit
    }
}

#[derive(Debug, Clone)]
pub struct Identity {
    qubit: usize,
}

impl Identity {
    pub fn new(qubit: usize) -> Self {
        Self { qubit }
    }
}

impl QuantumGate for Identity {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        array![
            [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
            [Complex::new(0.0, 0.0), Complex::new(1.0, 0.0)]
        ]
    }

    fn name(&self) -> String {
        String::from("I")
    }
}

impl SingleQubit for Identity {
    fn target_qubit(&self) -> usize {
        self.qubit
    }
}
