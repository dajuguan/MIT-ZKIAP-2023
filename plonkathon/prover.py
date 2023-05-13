from compiler.program import Program, CommonPreprocessedInput
from utils import *
from setup import *
from typing import Optional
from dataclasses import dataclass
from transcript import Transcript, Message1, Message2, Message3, Message4, Message5
from poly import Polynomial, Basis


@dataclass
class Proof:
    msg_1: Message1
    msg_2: Message2
    msg_3: Message3
    msg_4: Message4
    msg_5: Message5

    def flatten(self):
        proof = {}
        proof["a_1"] = self.msg_1.a_1
        proof["b_1"] = self.msg_1.b_1
        proof["c_1"] = self.msg_1.c_1
        proof["z_1"] = self.msg_2.z_1
        proof["t_lo_1"] = self.msg_3.t_lo_1
        proof["t_mid_1"] = self.msg_3.t_mid_1
        proof["t_hi_1"] = self.msg_3.t_hi_1
        proof["a_eval"] = self.msg_4.a_eval
        proof["b_eval"] = self.msg_4.b_eval
        proof["c_eval"] = self.msg_4.c_eval
        proof["s1_eval"] = self.msg_4.s1_eval
        proof["s2_eval"] = self.msg_4.s2_eval
        proof["z_shifted_eval"] = self.msg_4.z_shifted_eval
        proof["W_z_1"] = self.msg_5.W_z_1
        proof["W_zw_1"] = self.msg_5.W_zw_1
        return proof


@dataclass
class Prover:
    group_order: int
    setup: Setup
    program: Program
    pk: CommonPreprocessedInput

    def __init__(self, setup: Setup, program: Program):
        self.group_order = program.group_order
        self.setup = setup
        self.program = program
        self.pk = program.common_preprocessed_input()

    def prove(self, witness: dict[Optional[str], int]) -> Proof:
        # Initialise Fiat-Shamir transcript
        transcript = Transcript(b"plonk")

        # Collect fixed and public information
        # FIXME: Hash pk and PI into transcript
        public_vars = self.program.get_public_assignments()
        PI = Polynomial(
            [Scalar(-witness[v]) for v in public_vars]
            + [Scalar(0) for _ in range(self.group_order - len(public_vars))],
            Basis.LAGRANGE,
        )
        self.PI = PI   # public input polynomial

        # Round 1
        msg_1 = self.round_1(witness)
        self.beta, self.gamma = transcript.round_1(msg_1)

        # Round 2
        msg_2 = self.round_2()
        self.alpha, self.fft_cofactor = transcript.round_2(msg_2)

        # Round 3
        msg_3 = self.round_3()
        self.zeta = transcript.round_3(msg_3)

        # Round 4
        msg_4 = self.round_4()
        self.v = transcript.round_4(msg_4)

        # Round 5
        msg_5 = self.round_5()

        return Proof(msg_1, msg_2, msg_3, msg_4, msg_5)

    def round_1(
        self,
        witness: dict[Optional[str], int],
    ) -> Message1:
        program = self.program
        setup = self.setup
        group_order = self.group_order   #门约束的等式个数最接近的2^K（at least 8）
        print("self.pk.QL--------->",len(self.pk.QL.values))
        if None not in witness:
            witness[None] = 0

        # Compute wire assignments for A, B, C, corresponding:
        # - A_values: witness[program.wires()[i].L]
        # - B_values: witness[program.wires()[i].R]
        # - C_values: witness[program.wires()[i].O]
        A_values, B_values, C_values = [], [], []
        for i in range(group_order):
            if i < len(program.wires()):
                A_values.append(Scalar(1)*witness[program.wires()[i].L])
                B_values.append(Scalar(1)*witness[program.wires()[i].R])
                C_values.append(Scalar(1)*witness[program.wires()[i].O])
            else:
                A_values.append(Scalar(0))
                B_values.append(Scalar(0))
                C_values.append(Scalar(0))

        # Construct A, B, C Lagrange interpolation polynomials for
        # A_values, B_values, C_values
        self.A = Polynomial(A_values, Basis.LAGRANGE)
        self.B = Polynomial(B_values, Basis.LAGRANGE)
        self.C = Polynomial(C_values, Basis.LAGRANGE)

        # Compute a_1, b_1, c_1 commitments to A, B, C polynomials
        a_1 = setup.commit(self.A)
        b_1 = setup.commit(self.B)
        c_1 = setup.commit(self.C)

        # Sanity check that witness fulfils gate constraints
        assert (
            self.A * self.pk.QL
            + self.B * self.pk.QR
            + self.A * self.B * self.pk.QM
            + self.C * self.pk.QO
            + self.PI
            + self.pk.QC
            == Polynomial([Scalar(0)] * group_order, Basis.LAGRANGE)
        )
        print("Round 1:gate constraints fullfilled")
        # Return a_1, b_1, c_1
        return Message1(a_1, b_1, c_1)

    def round_2(self) -> Message2:
        group_order = self.group_order
        setup = self.setup
        roots_of_unity = Scalar.roots_of_unity(group_order)
        # Using A, B, C, values, and pk.S1, pk.S2, pk.S3, compute
        # Z_values for permutation grand product polynomial Z
        #
        # Note the convenience function:
        #       self.rlc(val1, val2) = val_1 + self.beta * val_2 + gamma
        Z_values = [Scalar(1)]
        for i in range(group_order):  # must be order+1 to calculate Z_n
            # f_i =  self.rlc(self.A.values[i-1], roots_of_unity[i-1])* self.rlc(self.B.values[i-1], 2 * roots_of_unity[i-1])* self.rlc(self.C.values[i-1], 3 * roots_of_unity[i-1])
            # g_i =  self.rlc(self.A.values[i-1], self.pk.S1.values[i-1]) * self.rlc(self.B.values[i-1], self.pk.S2.values[i-1]) * self.rlc(self.C.values[i-1], self.pk.S3.values[i-1])
            f_i =  self.rlc(self.A.values[i], roots_of_unity[i])* self.rlc(self.B.values[i], 2 * roots_of_unity[i])* self.rlc(self.C.values[i], 3 * roots_of_unity[i])
            g_i =  self.rlc(self.A.values[i], self.pk.S1.values[i]) * self.rlc(self.B.values[i], self.pk.S2.values[i]) * self.rlc(self.C.values[i], self.pk.S3.values[i])
            z_i_plus_1 = f_i*Z_values[-1]/g_i
            Z_values.append(Scalar(z_i_plus_1))
            print("z_i_plus_1------>",z_i_plus_1)

        # Check that the last term Z_n = 1
        assert Z_values.pop() == 1

        # Sanity-check that Z was computed correctly
        for i in range(group_order):
            assert (
                self.rlc(self.A.values[i], roots_of_unity[i])
                * self.rlc(self.B.values[i], 2 * roots_of_unity[i])
                * self.rlc(self.C.values[i], 3 * roots_of_unity[i])
            ) * Z_values[i] - (
                self.rlc(self.A.values[i], self.pk.S1.values[i])
                * self.rlc(self.B.values[i], self.pk.S2.values[i])
                * self.rlc(self.C.values[i], self.pk.S3.values[i])
            ) * Z_values[
                (i + 1) % group_order
            ] == 0
        print("Round 2:copy constraints fullfilled")

        # Construct Z, Lagrange interpolation polynomial for Z_values
        # Cpmpute z_1 commitment to Z polynomial
        self.Z = Polynomial(Z_values,Basis.LAGRANGE)
        z_1 = setup.commit(self.Z)

        # Return z_1
        return Message2(z_1)

    def round_3(self) -> Message3:
        group_order = self.group_order
        setup = self.setup
        fft_cofactor = self.fft_cofactor

        # Compute the quotient polynomial

        # List of roots of unity at 4x fineness, i.e. the powers of µ
        # where µ^(4n) = 1
        quarter_roots = Scalar.roots_of_unity(group_order * 4)
        # Using self.fft_expand, move A, B, C into coset extended Lagrange basis
        A = self.fft_expand(self.A)
        B = self.fft_expand(self.B)
        C = self.fft_expand(self.C)

        # Expand public inputs polynomial PI into coset extended Lagrange
        PI = self.fft_expand(self.PI)

        # Expand selector polynomials pk.QL, pk.QR, pk.QM, pk.QO, pk.QC
        # into the coset extended Lagrange basis
        QL = self.fft_expand(self.pk.QL)
        QR = self.fft_expand(self.pk.QR)
        QM = self.fft_expand(self.pk.QM)
        QO = self.fft_expand(self.pk.QO)
        QC = self.fft_expand(self.pk.QC)

        # Expand permutation grand product polynomial Z into coset extended
        # Lagrange basis
        Z = self.fft_expand(self.Z)

        # Expand shifted Z(ω) into coset extended Lagrange basis
        # Zw = self.fft_expand(self.Zw)  #Z(wx)
        Zw = Z.shift(4)
        # Expand permutation polynomials pk.S1, pk.S2, pk.S3 into coset
        # extended Lagrange basis
        S1 = self.fft_expand(self.pk.S1)
        S2 = self.fft_expand(self.pk.S2)
        S3 = self.fft_expand(self.pk.S3)

        # Compute Z_H = X^N - 1, also in evaluation form in the coset
        Z_H = Polynomial(
            [
                ((Scalar(r) * fft_cofactor) ** group_order -1)  for r in quarter_roots
            ],
            Basis.LAGRANGE
        )

        self.Z_H = Z_H
        quarter_roots = Polynomial(
            Scalar.roots_of_unity(group_order * 4), Basis.LAGRANGE
        )
        # Compute L0, the Lagrange basis polynomial that evaluates to 1 at x = 1 = ω^0
        # and 0 at other roots of unity

        # Expand L0 into the coset extended Lagrange basis
        L0_big = self.fft_expand(
            Polynomial([Scalar(1)] + [Scalar(0)] * (group_order - 1), Basis.LAGRANGE)
        )
        self.L0_big = L0_big

        # Compute the quotient polynomial (called T(x) in the paper)
        # It is only possible to construct this polynomial if the following
        # equations are true at all roots of unity {1, w ... w^(n-1)}:
        # 1. All gates are correct:
        #    A * QL + B * QR + A * B * QM + C * QO + PI + QC = 0
        #
        # 2. The permutation accumulator is valid:
        #    Z(wx) = Z(x) * (rlc of A, X, 1) * (rlc of B, 2X, 1) *
        #                   (rlc of C, 3X, 1) / (rlc of A, S1, 1) /
        #                   (rlc of B, S2, 1) / (rlc of C, S3, 1)
        #    rlc = random linear combination: term_1 + beta * term2 + gamma * term3
        #
        # 3. The permutation accumulator equals 1 at the start point
        #    (Z - 1) * L0 = 0
        #    L0 = Lagrange polynomial, equal at all roots of unity except 1

        QUOT_big = ( 
            (A*B*QM + A*QL + B*QR + C*QO + PI + QC) +
            (self.rlc(A,quarter_roots * fft_cofactor) * self.rlc(B,quarter_roots * (fft_cofactor * 2)) * self.rlc(C,quarter_roots * (fft_cofactor * 3)) * Z * self.alpha) -
            (self.rlc(A,S1) * self.rlc(B,S2) * self.rlc(C,S3) * Zw * self.alpha ) +
            (Z - Scalar(1)) * L0_big * self.alpha**2 
        )
        QUOT_big = QUOT_big / Z_H
        self.QUOT_big = ( 
            (A*B*QM + A*QL + B*QR + C*QO + PI + QC)
            # (self.rlc(A,quarter_roots * fft_cofactor) * self.rlc(B,quarter_roots * (fft_cofactor * 2)) * self.rlc(C,quarter_roots * (fft_cofactor * 3)) * Z * self.alpha) -
            # (self.rlc(A,S1) * self.rlc(B,S2) * self.rlc(C,S3) * Zw * self.alpha ) +
            # (Z - Scalar(1)) * L0_big * self.alpha**2 
        )/ Z_H

        # Sanity check: QUOT has degree < 3n
        assert (
            self.expanded_evals_to_coeffs(QUOT_big).values[-group_order:]
            == [0] * group_order
        )
        print("Generated the quotient polynomial------------->")

        # Split up T into T1, T2 and T3 (needed because T has degree 3n - 4, so is
        # too big for the trusted setup)
        coeffs = self.expanded_evals_to_coeffs(QUOT_big).values

        T1 = Polynomial(coeffs[:group_order],Basis.MONOMIAL).fft()
        T2 = Polynomial(coeffs[group_order:2*group_order],Basis.MONOMIAL).fft()
        T3 = Polynomial(coeffs[2*group_order: group_order * 3],Basis.MONOMIAL).fft()
        self.t_lo = T1
        self.t_mid = T2
        self.t_hi = T3

        # Sanity check that we've computed T1, T2, T3 correctly
        assert (
            T1.barycentric_eval(fft_cofactor)
            + T2.barycentric_eval(fft_cofactor) * fft_cofactor**group_order
            + T3.barycentric_eval(fft_cofactor) * fft_cofactor ** (group_order * 2)
        ) == QUOT_big.values[0]

        print("Generated T1, T2, T3 polynomials")
        print("QUOT_big.values[0]----->",QUOT_big.values[0])

        # Compute commitments t_lo_1, t_mid_1, t_hi_1 to T1, T2, T3 polynomials
        t_lo_1 = setup.commit(T1)
        t_mid_1 = setup.commit(T2)
        t_hi_1 = setup.commit(T3)
        # print("t_lo_1-------->",t_lo_1)
        # print("t_lo_2-------->",t_mid_1)
        # print("t_hi_1-------->",t_hi_1)
        # Return t_lo_1, t_mid_1, t_hi_1
        return Message3(t_lo_1, t_mid_1, t_hi_1)

    def round_4(self) -> Message4:
        # Compute evaluations to be used in constructing the linearization polynomial.
        zeta = self.zeta
        group_order = self.group_order
        root_of_unity = Scalar.root_of_unity(group_order)
        # Compute a_eval = A(zeta)
        # Compute b_eval = B(zeta)
        # Compute c_eval = C(zeta)
        # Compute s1_eval = pk.S1(zeta)
        # Compute s2_eval = pk.S2(zeta)
        # Compute z_shifted_eval = Z(zeta * ω)
        a_eval = self.A.barycentric_eval(zeta)
        b_eval = self.B.barycentric_eval(zeta)
        c_eval = self.C.barycentric_eval(zeta)
        s1_eval = self.pk.S1.barycentric_eval(zeta)
        s2_eval = self.pk.S2.barycentric_eval(zeta)
        z_shifted_eval = self.Z.barycentric_eval(zeta * root_of_unity)
        # print("a_eval------->",a_eval)
        
        self.a_eval = a_eval
        self.b_eval = b_eval
        self.c_eval = c_eval
        self.s1_eval = s1_eval
        self.s2_eval = s2_eval
        self.z_shifted_eval = z_shifted_eval

        # Return a_eval, b_eval, c_eval, s1_eval, s2_eval, z_shifted_eval
        return Message4(a_eval, b_eval, c_eval, s1_eval, s2_eval, z_shifted_eval)

    def round_5(self) -> Message5:
        group_order = self.group_order
        setup = self.setup
        zeta = self.zeta

        # Evaluate the Lagrange basis polynomial L0 at zeta
        L0_ev = Polynomial(
            [Scalar(1)] + [Scalar(0)] * (group_order - 1), Basis.LAGRANGE
        ).barycentric_eval(zeta)
        print("L0_ev----->", L0_ev)
        print("L0_ev----->", self.L0_big.barycentric_eval(zeta/self.fft_cofactor))  #经过了expand，所以计算的时候也需要烦着expand回来
        assert(L0_ev == self.L0_big.barycentric_eval(zeta/self.fft_cofactor))

        # Evaluate the vanishing polynomial Z_H(X) = X^n - 1 at zeta
        # Z_H_ev = self.Z_H.barycentric_eval(zeta)
        Z_H_ev = zeta**(group_order) - 1

        # Move T1, T2, T3 into the coset extended Lagrange basis
        # Move pk.QL, pk.QR, pk.QM, pk.QO, pk.QC into the coset extended Lagrange basis
        # Move Z into the coset extended Lagrange basis
        # Move pk.S3 into the coset extended Lagrange basis
        t_lo = self.fft_expand(self.t_lo)   # no need to do fft_expand, because the degree will not exceed n
        t_mid = self.fft_expand(self.t_mid)
        t_hi = self.fft_expand(self.t_hi)
        QL = self.fft_expand(self.pk.QL)
        QR = self.fft_expand(self.pk.QR)
        QM = self.fft_expand(self.pk.QM)
        QO = self.fft_expand(self.pk.QO)
        QC = self.fft_expand(self.pk.QC)
        Z = self.fft_expand(self.Z)
        S3 = self.fft_expand(self.pk.S3)

        # Compute the "linearization polynomial" R. This is a clever way to avoid
        # needing to provide evaluations of _all_ the polynomials that we are
        # checking an equation betweeen: instead, we can "skip" the first
        # multiplicand in each term. The idea is that we construct a
        # polynomial which is constructed to equal 0 at Z only if the equations
        # that we are checking are correct, and which the verifier can reconstruct
        # the KZG commitment to, and we provide proofs to verify that it actually
        # equals 0 at Z
        #
        # In order for the verifier to be able to reconstruct the commitment to R,
        # it has to be "linear" in the proof items, hence why we can only use each
        # proof item once; any further multiplicands in each term need to be
        # replaced with their evaluations at Z, which do still need to be provided

        # Commit to R
        PI_ev = self.PI.barycentric_eval(zeta)
        print("len(self.PI.values)",len(self.PI.values))
        alpha = self.alpha
        s1_eval = self.s1_eval
        s2_eval = self.s2_eval
        z_shifted_eval = self.z_shifted_eval

        c_eval = Polynomial([self.c_eval] * group_order, Basis.LAGRANGE)  #由于只能是polynomial*Scalar因此需要将c_eval拓展

        R = (
            self.pk.QM*self.a_eval*self.b_eval + self.pk.QL* self.a_eval + self.pk.QR*self.b_eval + self.pk.QO*self.c_eval + PI_ev + self.pk.QC +
            (
                self.Z*self.rlc(self.a_eval, zeta)*self.rlc(self.b_eval, zeta*2)*self.rlc(self.c_eval, zeta*3)
                - self.rlc(c_eval, self.pk.S3)*self.rlc(self.a_eval, s1_eval)*self.rlc(self.b_eval, s2_eval)*z_shifted_eval
            )*alpha + 
            (self.Z - Scalar(1))*L0_ev*(alpha**2)
            - (self.t_lo + self.t_mid*(zeta**group_order) + self.t_hi*(zeta**(2*group_order)))*Z_H_ev
        )
        
        # the bellow is the same as the above code without expanding
        # c_eval = Polynomial([self.c_eval] * group_order * 4, Basis.LAGRANGE)  #由于只能是polynomial*Scalar因此需要将c_eval拓展
        
        # R_big = (
        #     QM*self.a_eval*self.b_eval + QL* self.a_eval + QR*self.b_eval + QO*self.c_eval + PI_ev + QC +
        #     (
        #         Z*self.rlc(self.a_eval, zeta)*self.rlc(self.b_eval, zeta*2)*self.rlc(self.c_eval, zeta*3)
        #         - self.rlc(c_eval, S3)*self.rlc(self.a_eval, s1_eval)*self.rlc(self.b_eval, s2_eval)*z_shifted_eval
        #     )*alpha + 
        #     (Z - Scalar(1))*L0_ev*(alpha**2)
        #     - (t_lo + t_mid*(zeta**group_order) + t_hi*(zeta**(2*group_order)))*Z_H_ev
        # )

        # R_coeffs = self.expanded_evals_to_coeffs(R_big).values
        # assert R_coeffs[group_order:] == [0] * (group_order * 3)
        # R = Polynomial(R_coeffs[:group_order], Basis.MONOMIAL).fft()

        # Sanity-check R
        assert R.barycentric_eval(zeta) == 0

        print("Generated linearization polynomial R")

        # Generate proof that W(z) = 0 and that the provided evaluations of
        # A, B, C, S1, S2 are correct

        # Move A, B, C into the coset extended Lagrange basis
        # Move pk.S1, pk.S2 into the coset extended Lagrange basis

        # In the COSET EXTENDED LAGRANGE BASIS,
        # Construct W_Z = (
        #     R
        #   + v * (A - a_eval)
        #   + v**2 * (B - b_eval)
        #   + v**3 * (C - c_eval)
        #   + v**4 * (S1 - s1_eval)
        #   + v**5 * (S2 - s2_eval)
        # ) / (X - zeta)
        v = self.v
        W_Z = (
            R
            + (self.A - self.a_eval)*v 
            + (self.B - self.b_eval)*v**2
            + (self.C - self.c_eval)*v**3
            + (self.pk.S1 - s1_eval)*v**4
            + (self.pk.S2 - s2_eval)*v**5
        ) / Polynomial([-zeta, Scalar(1)] + [Scalar(0)]*(group_order - 2), Basis.MONOMIAL).fft()
        # W_z_coeffs = self.expanded_evals_to_coeffs(W_Z).values
        # Check that degree of W_z is not greater than n
        # assert W_z_coeffs[group_order:] == [0] * (group_order * 3)

        # Compute W_z_1 commitment to W_z
        W_z_1 = setup.commit(W_Z)

        # Generate proof that the provided evaluation of Z(z*w) is correct. This
        # awkwardly different term is needed because the permutation accumulator
        # polynomial Z is the one place where we have to check between adjacent
        # coordinates, and not just within one coordinate.
        # In other words: Compute W_zw = (Z - z_shifted_eval) / (X - zeta * ω)
        ω = Scalar.root_of_unity(group_order)
        W_zw = (self.Z - z_shifted_eval) / Polynomial([-zeta*ω, Scalar(1)] + [Scalar(0)]*(group_order - 2), Basis.MONOMIAL).fft()
        W_zw_coeffs = self.expanded_evals_to_coeffs(W_zw).values
        # Check that degree of W_z is not greater than n
        # assert W_zw_coeffs[group_order:] == [0] * (group_order * 3)

        # Compute W_zw_1 commitment to W_z
        W_zw_1 = setup.commit(W_zw)
        print("Generated final quotient witness polynomials")

        # Return W_z_1, W_zw_1
        return Message5(W_z_1, W_zw_1)

    def fft_expand(self, x: Polynomial):
        return x.to_coset_extended_lagrange(self.fft_cofactor)

    def expanded_evals_to_coeffs(self, x: Polynomial):
        return x.coset_extended_lagrange_to_coeffs(self.fft_cofactor)

    def rlc(self, term_1, term_2):
        return term_1 + term_2 * self.beta + self.gamma
