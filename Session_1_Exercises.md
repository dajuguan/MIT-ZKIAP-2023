# ZKP for 3-coloring Demo
## Exercise 1
Q: Currently, you can only select adjacent pairs of nodes to check. Would the proof still be zero knowledge if you could pick arbitrary pairs of nodes to check?

A: It's not. Because it may violate the Soundness properties of ZKP, if the prover provide the wrong answer, he may not be caught.

## Exercise 2
Q: The equation currently being used for confidence is 1-(1/E)^n, where E is the number of edges in the graph, and n is the number of trials run. Is this the correct equation? Why is there no prior?

A: Yes. Because the prover can change colors between different rounds of the game, so the event of picking an edge each tiem is independent.

# Optional - ZKP for DLOG

```python ./Session_1_code.py```

# zkmessage.xyz
## Q1
Q: Explain why you need to generate and save a “secret” value.

A: We will use the "secret" value to claim a ZKMessage public key and generate group signatures.

## Q2
Q: Write out a plain-English explanation of what statement is being proven in ZK.

A: I know some private key that represents the author of the post belongs to a group of people listed in the group-signed message. Here is a signature that I know such a private key, without telling you what the acture private key is.

## Q3
Q: Log into the same zkmessage account, from a different browser or computer. Explain why zkmessage can’t just use a simple “username/password” system like most social apps.

A: In that way, we have to give the control of our account to the central server, thus they delete our account or messgages. But if we use the private key to login, we can control our account. If the server  censor our messages or signatures, we can swap out the server for a decentralized data layer.