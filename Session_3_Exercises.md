# Quadratic nonresidue
## (a) Completeness:
Q:  if QR(m,x) = 0 and both parties behave according to the protocol, then the Verifier always accepts.

A: There are two cases to check.
1. if b = 0, we should check that $ s^2x$ is not a quadratic residue mod m. Suppose $QR(m,s^2x)=1$ , i.e., there exists an integer $r$ such that $r^2 \equiv s^2x \:(mod\:m)$. Because $s$ is relatively prime to $m$, there must exist an integer $z$ such that $(r^2 - s^2z^2) \equiv 0 \:(mod\:m)$. So $z^2 \equiv x \:(mod\:m)$, it violates the assumption that $QR(m,x)=0$. Thus, if b=0, then $QR(m,y)=0$.

2. if b = 1, it's easy to see there exists an integer $s$ such that  $QR(m,y) = 1$.

Hence, the completeness is ensured if $QR(m,s)=0$.

## (b) Soundness: 
Q: if QR(m,x) = 1, then no matter what the Prover does (the Prover does not have to follow the protocol), the Verifier rejects with probability ≥ 1/2.

A: Base on question (a)'s answer, if $QR(m,x) = 1$, then $QR(m,y) = 1$ no matter b=0 or 1. Then, there are two cases.
1. If the Prover follows the protocol, when b=0, he sends 1 to Verifier, the Verifier rejects because 1 is not equal to 0; when b=1, the Verifier will accept in the same way. Hence, the Verifier rejects with probability  1/2.

2. If the Prover totally doesn't follow the protocol, he will send 0 to Verifier no matter b=0 or 1. So the Verifier will accept when b=0, and reject when b=1. Hence, the Verifier rejects with probability  1/2.

Hence, no matter what the Prover does,the Verifier rejects with probability ≥ 1/2.

# Quadratic residue
## (a) Completeness:
Q:  If both parties behaves according to the protocol, then the Verifier always accepts.

A: The Verifier will always accepts because
$$
y \equiv
   \begin{cases}
     u^2x = xt^2, \quad if\: b=0\\
     s^2t^2=xt^2, \quad if \:  b=1
   \end{cases}
$$

## (b) Soundness:
Q: If QR(m, x) = 0 (in particular, the Prover does not know any valid s), then no matter what the Prover does, the Verifer rejects with probability ≥ 1/2.

A: When b=0, the Verifier accepts because $y \equiv xt^2 $. However, when b=1, $y\equiv s^2t^2\neq xt^2 \:(mod \: m)$, so the Verifier rejects.

## (c) Zero knowledge
Q: No matter what the Verifier does (the Verifier could be- have differently from the Protocol), the Verifier could have simulated the entire interaction on its own without interacting with the Prover and such that if the Verifier accepts, then the transcript is indistinguishable from the actual interaction.

A: If the Verifier chooses b=0, then he gets two messages: $(xt^2, t)$, which doesn't reveal anything about $s$. If the Verifier chooses b=1, then he gets two messages: $(xt^2, st)$ which looks like two random messages. In both cases, these two messages are completely independent of $s$, and hence reveals no additional information.

# Self-pairing implies failure of DDH

Because $e(g,g)^{(αβ)} = e(g^α,g^β)$, and calculating the function $e$ is efficient, we can just test $e(g^α,g^β) = e(g,g^y)$ to determine whether $αβg = y$.

[1] 张方国.从双线性对到多线性映射[J].密码学报,2016,3(03):211-228.[DOI:10.13868/j.cnki.jcr.000122](https://kns.cnki.net/kcms2/article/abstract?v=3uoqIhG8C44YLTlOAiTRKibYlV5Vjs7ijP0rjQD-AVm8oHBO0FTadjnbZbZqJIgs1Xow4kTBM7qIkSmHbzozfuIzhPf1mNsY&uniplatform=NZKPT).

# BLS signature aggregation

# Q1
Q: Check that the verification algorithm always accepts a correctly signed signature.

A: Based on the properties of Pairing, we can derive:
$$
e(g^0,σ) = e(g^0,αH(m)) = e(αg^0,H(m))= e(pk,H(m))
$$
Thus, the verification algorithm always accepts a correctly signed signature.

# Q2
Q: Can you specify some computational hardness assumptions?

A: Security of BLS relies on the standard **co-CDH assumption** in the bilinear group $(G_0, G_1)$. The assumption states that for all efficient algorithms $\Alpha$ the quantity 
$$
CDHadv[\Alpha, (G_0, G_1)]\: := Pr[\Alpha(g_1^\alpha, g_0^\beta, g_0^\alpha)=g_0^{\alpha\beta}]
$$
is negligible, where $\alpha, \beta 
\stackrel{R}{\leftarrow}Z_q$.

[Ref: BLS Multi-Signatures With Public-Key Aggregation](https://crypto.stanford.edu/~dabo/pubs/papers/BLSmultisig.html)
