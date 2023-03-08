## Q1
Q: The Setup phase of the KZG polynomial commitment scheme involves computing commitments to powers of a secret evaluation point τ . This is called the “trusted setup" and is often generated in a multi-party computation known as the “Powers of Tau" ceremony. One day, you find the value of τ on a slip of paper. How can you use it to make a fake KZG opening proof?

A: In the KZG open phase, the verifier needs to check:
$$
e(π, [τ − x]_2) = e(C − [y]_1, G_2) 
$$
$$
{f(τ )-y \over τ  -x}(τ − x)= f(τ) -y
$$
so, it's easy to see we can construct arbitrary $y'=f(x')$ to make a fake 
$$
τ'= {f(τ )-y' \over τ  -x}
$$

## Q2
Q: Construct a vector commitment scheme from the KZG polynomial commitment scheme.
(Hint: For a vector m = (m1, . . . , mq), is there an “interpolation polynomial" I(X) such
that I(i) = m[i]?)

A:
We can construct a Lagrange interpolation polynomial I(X) such that:
$$
I(X) = \sum_i m_i L_i(X) = \left\{
\begin{aligned}
m_i \:\: if X=i \\
0 \:\: else
\end{aligned}
\right.
$$

- KZG.Setup: the same as KZG commitment scheme
- KZG.Commit $(ck; I(X)) ->C$: for $C=I(\alpha)_1$
- KZG.Openn $(srs, C, m_i, i) → {0, 1}$: To “open" the commitment at evaluation point x
to a claimed value y.
    - the procedure is the same as the KZG commitment scheme


## Q3
Q:The KZG polynomial commitment scheme makes an opening proof π for the relation p(x) =
y. Can you extend the scheme to produce a multiproof π, that convinces us of p(xi) = yi
for a list of points and evaluations (xi
, yi)? (Hint: assume that you have an interpolation
polynomial I(X) such that I(xi) = yi.

A: Construct an interpolation polynomial:
$$
I(X) = \sum_i y_i L_i(X) = \left\{
\begin{aligned}
y_i \:\: if X=i \\
0 \:\: else
\end{aligned}
\right.
$$

To extend the KZG proof to multiple points and evaluations ${(x_i, y_i)}$, we have to ensure every $x_i$ is the root of the numerator of $π_{multiproof}$. So, we can set:

$$
π_{multiproof} := {f(X) - I(X) \over Π(X-x_i)}
$$

And the verifier should check:
$$
e(π_{multiproof}, [Π(X-x_i)]_2) = e(C − I(X), G_2) 
$$
