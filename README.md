# Pollard's Lambda Algorithm for SECP256K1 Curve (ρλ)

![C++](https://img.shields.io/badge/Language-C++-blue)
![Linux](https://img.shields.io/badge/Platform-Linux-white)

## Fast Links:

[Main Index](https://github.com/lucaselblanc/pollardslambda/tree/main?tab=readme-ov-file#pollards-lambda-algorithm-for-secp256k1-curve-%CF%81%CE%BB)

[Description](https://github.com/lucaselblanc/pollardslambda/tree/main?tab=readme-ov-file#description)

[Benchmark](https://github.com/lucaselblanc/pollardslambda/tree/main#benchmark-tpu-v5e-8)

[Average K Factor](https://github.com/lucaselblanc/pollardslambda/tree/main#average-k-factor)

[Technical Features](https://github.com/lucaselblanc/pollardslambda/tree/main#technical-features)

[Distinguished Points](https://github.com/lucaselblanc/pollardslambda/tree/main#distinguished-points-dp)

[Delay Of Distinguished Points](https://github.com/lucaselblanc/pollardslambda/tree/main#delay-of-distinguished-points)

[Algorithm Complexity](https://github.com/lucaselblanc/pollardslambda/tree/main#algorithm-complexity)

[Negation Map](https://github.com/lucaselblanc/pollardslambda/tree/main#negation-map)

[Academic Paper References](https://github.com/lucaselblanc/pollardslambda/tree/main#academic-references)

[Prerequisites](https://github.com/lucaselblanc/pollardslambda/tree/main#prerequisites)

[Installation](https://github.com/lucaselblanc/pollardslambda/tree/main#installation)

[Commands](https://github.com/lucaselblanc/pollardslambda/tree/main#commands)

[External Libraries Used](https://github.com/lucaselblanc/pollardslambda/tree/main#external-libraries-used)

## Notice

 This project is still full CPU, migration to GPU is being analyzed.

## Description

 This repository contains a high-performance implementation of Pollard’s Lambda algorithm for solving the Elliptic Curve Discrete Logarithm Problem (ECDLP) on the secp256k1 curve.

#### Pollard's Lambda (ρλ)

 The algorithm uses parallel pseudo-random walks based on an R-adding
walk construction.

Each walker maintains:

$R = a \cdot G + b \cdot H$

`G` is the generator point.

`H` is the target public key.

`a` and `b` are scalar coefficients used for recovery after collision.

The transition function uses a precomputed step table and MurmurHash3
avalanche mixing.

When two walkers reach the same point with different coefficients, a
modular equation allows recovery through modular inversion.

## Benchmark CPU v5e-8 224 cores

<details>
<summary><strong>  Note:</strong></summary>

<sub><i>Performance metrics are governed by probabilistic Monte Carlo variables. Consequently, execution cycles will demonstrate natural stochastic fluctuation across independent iterations, asymptotically converging upon the expected mean K-factor of the algorithm.</i></sub>

</details>

```
5 bits ≈ 00:00:00
10 bits ≈ 00:00:00
15 bits ≈ 00:00:00
20 bits ≈ 00:00:00
25 bits ≈ 00:00:00
30 bits ≈ 00:00:00
35 bits ≈ 00:00:00
40 bits ≈ 00:00:00
45 bits ≈ 00:00:00
50 bits ≈ 00:00:00
55 bits ≈ 00:00:04
60 bits ≈ 00:00:07
65 bits ≈ 00:00:43
70 bits ≈ 00:01:09
```

## Average k-Factor

The formal equation by van Oorschot and Wiener (1999) on paper: [P. C. Van Oorschot & M. J. Wiener (1999) - Parallel Collision Search with Cryptanalytic Applications](https://people.scs.carleton.ca/~paulv/papers/JoC97.pdf) for the expected total time of the Lambda method operating on a restricted interval is:

$T_\lambda = \left( \frac{2 \sqrt{b}}{m} + \frac{1}{\theta} \right) t$

The expected constant for the factor k is:

$E[\text{ops}] \approx 1.2533 \sqrt{W}$

This value represents the theoretical efficiency boundary predicted for parallel collision search when elliptic curve symmetry reduction is applied through the Negation Map.

The implementation achieves an empirical average k-factor:

<details>
<summary><strong>  Note:</strong></summary>

<sub><i>The project is under active development. The reported average k-factor reflects the algorithm's current state and should not be considered the final value for implementation, and the actual average k-factor of this implementation will never be achieved in smaller intervals because in the context of Monte Carlo-based collision search algorithms for ECDLP, the average K-factor exhibits deterministic asymptotic stabilization as the search space (key range) expands. This performance convergence is intrinsically driven by the progressive vanishing of geometric boundary penalties. In exponentially larger intervals, spatial density dilutes, causing the probability of walker ejections or 'deaths' resulting from edge collisions to approach zero, thereby eradicating fruitless cycles. Consequently, the system's optimal mean step size multiplier ($m$) organically relaxes toward the continuous theoretical limit of 3.0, corroborating established cryptographic literature on random walks. This thermodynamic behavior is strictly governed by Markov chain dynamics, wherein the initial boundary-induced clustering penalty undergoes a probabilistic exponential decay, ultimately converging into its asymptotic stationary distribution.</i></sub>

</details>

```cpp
//8 thousand samples in the 50-bit range.
```

$average k ≈ 1.230$

$median k ≈ 1.157$

This value was obtained through thousands of independent benchmark samples, measuring the average number of operations required to solve interval-restricted searches.

The observed k-factor is the result of the combined effect of multiple engineering optimizations:

— Negation Map optimization.

— Parallel collision search architecture.

— Optimized random walk distribution.

— Distinguished Points collision detection.

— Cache-aware precomputed step windows.

— Efficient walker synchronization.

— Batch Jacobian-to-Affine conversion.

Shows that the implementation operates below the theoretical average bound expected for the generic optimized Pollard's Lambda model, approaching practical optimal performance under real execution conditions.

Lower-than-average k-factor values may still occur due to the statistical nature of random walks, where some searches benefit from favorable trajectory convergence and collision timing.

## Implementation Highlights

— Multi-threaded walkers.

— Cache-aware tables.

— Batch inversion.

— Snapshot recovery.

— Distinguished Point system.

— Negation Map optimization.

## Technical Features

#### Batch Jacobian-to-Affine (Montgomery Trick)

 Inversions in finite fields are computationally expensive. Both versions utilize Batch Inversion, processing multiple walkers simultaneously. This allows the algorithm to perform only one modular inversion per batch, converting Jacobian coordinates to Affine at a fraction of the usual cost, only the x coordinate is calculated to maintain efficiency.

#### Snapoints / Resumable Search State

 The entire search state: walker positions, scalar coefficients, and the
distinguished-point table, can be saved to disk and restored exactly where it left off. Saves are written atomically via a PID-namespaced temporary file promoted by a single `rename(2)`, with `fsync` on both the file and its parent directory before commit, guaranteeing full recovery even after a hard power loss. Every snapshot carries an integrity checksum and is validated against the current run parameters before any state is touched, so neither corruption nor a configuration mismatch can produce a silent incorrect resume.

#### Pre-Computed Points ```windowSize``` in L1/L2/L3 Caches

 The step window is dynamically adjusted according to the processor's cache behavior. While small searches benefit from larger entropy tables, large searches prioritize lower latency cache access, making the ultimate objective to perfectly balance step entropy and memory locality.

#### Distinguished Points (DP)

 The Distinguished Points strategy is a memory-saving filter. Instead of storing every step of the walk (which would crash your RAM), the algorithm only saves points that satisfy a specific condition: the first d bits of the x coordinate must be zero. When two walkers hit the same DP, a collision is found and the private key is recovered.
 
​The Trade-off:

More DP bits:

— Lower RAM usage.

— Slower detection.

Fewer DP bits:

— Faster detection.

— Higher RAM usage.

#### Delay Of Distinguished Points

 When a walker begins traversing a path already explored by another walker, a collision will be delayed if the distinguished points filter condition is not met for both walkers. The delay will be overcome after the distinct points are recorded in the dp table. The higher the dp value, the greater the delay for a collision to be detected and recorded by the hashmap, this inherent delay directly affects the algorithm's empirical k-factor. While the actual path convergence is dictated by the birthday paradox, the delayed detection forces walkers to perform additional operations before the collision is logged. Since the k-factor measures the total steps taken against the theoretical expectation, this DP overhead naturally inflates the final result.. To mitigate this, it would be necessary to disable the dp filter, but this would cause excessive RAM usage and would not be worth the effort, and would ruin the performance. This delay is a necessary evil when using distinct points.

Theoretical Calculus:

```
unsigned long long RAM_BYTES = (unsigned long long)sysconf(_SC_PAGESIZE) * (unsigned long long)sysconf(_SC_AVPHYS_PAGES);

int dp = std::max(1, std::max((int)std::ceil((key_range / 2.0) - std::log2((double)RAM_BYTES / 128.0 /*128 bytes*/)), key_range >> 2));
```

Simple Abstraction:

```
int dp = std::max<int>(1, std::min<int>(key_range >> 2, static_cast<int>(sizeof(int32_t) * CHAR_BIT)));
```

## Algorithm Complexity

 The expected time complexity of Pollard's Lambda algorithm for elliptic curves is $O(\sqrt{w})$, where $w$ is the width of the search interval. In this implementation, the probability distribution in the steps is restricted to $O\left(\sqrt{\frac{K}{2}}\right)$ by strictly bounding the random walk step sizes, which keeps the probabilistic window confined to the target range. Given secp256k1 and the inclusion of the Negation Map optimization, this translates to approximately $O\left(\sqrt{\frac{\text{range}}{2}}\right)$, as predicted by the birthday paradox for bounded random walks over a finite group.

## Walker Architecture

Each walker stores its current elliptic curve position and scalar state. Multiple independent walkers explore different trajectories simultaneously. This design enables large-scale CPU parallelization.

## Negation Map

The implementation applies elliptic curve negation symmetry (Equivalence Class Size 2).

$P = (x, y)$

$-P = (x, -y \bmod p)$

Equivalent points can be treated as the same state during the search.

The expected improvement changes the effective complexity to
approximately:

$O\left(\sqrt{\frac{\text{range}}{2}}\right)$

 This implementation enforces strict geometric bounds (2S for type 0 walks and 3S for type 2 walks) to prevent long tails and wasted CPU cycles. By cutting off extreme statistical bad luck scenarios at the 3S mark, the engine consistently achieves the expected theoretical average of k ≈ 1.2533.

 ​Due to the nature of the search, it may be possible to obtain solutions with a k-factor < 1.2533. These represent scenarios of extreme statistical luck solving the ECDLP within a distance of just 1S jump. This rare "sniper" event requires three specific conditions to align: the type 2 walker drops extremely close to the private key, immediately merges with a type 1 trail, and quickly triggers a Distinguished Point (DP). While rarer than standard 2S or 3S jumps convergences, the architecture fully capitalizes on these optimal drops when they occur.

## Academic References

[J. M. Pollard - Monte Carlo methods for index computation (mod p) (1978)](https://www.ams.org/journals/mcom/1978-32-143/S0025-5718-1978-0491431-9/S0025-5718-1978-0491431-9.pdf)

[Richard P. Brent - An improved Monte Carlo factorization algorithm (1980)](https://maths-people.anu.edu.au/~brent/pd/rpb051i.pdf)

[Peter L. Montgomery - Speeding the Pollard and Elliptic Curve Methods of Factorization (1987)](https://www.ams.org/journals/mcom/1987-48-177/S0025-5718-1987-0866113-7/S0025-5718-1987-0866113-7.pdf)

[P. C. Van Oorschot & M. J. Wiener - Parallel Collision Search with Cryptanalytic Applications (1999)](https://people.scs.carleton.ca/~paulv/papers/JoC97.pdf)

[Gaudry & Schost - A low-memory parallel version of Matsuo, Chao and Tsujii's algorithm (2004)](https://cs.uwaterloo.ca/~eschost/publications/ants.pdf)

[Galbraith & Ruprai - Using Equivalence Classes to Accelerate Solving
the Discrete Logarithm Problem in a Short
Interval (2010)](https://eprint.iacr.org/2010/615.pdf)
                                      
[SECG - SEC 2: Recommended Elliptic Curve Domain Parameters (secp256k1 specification) (2010)](https://www.secg.org/sec2-v2.pdf)

## Prerequisites

- g++
- boost/multiprecision/cpp_int.hpp
- libssl-dev

## Installation

1. Clone this repository:
    ```bash
    ~/$ git clone https://github.com/lucaselblanc/pollardslambda.git
    ```

2. Install the necessary libraries:
    ```bash
    sudo apt update
    sudo apt install build-essential -y
    sudo apt install boost-headers -y
    sudo apt install libssl-dev -y
    ```

3. Compile the project:
    ```bash
    ~/$cd pollardslambda
    ```

    ```bash
    ~/pollardslambda$ make
    ```

4. Run the program:
    ```bash
    ~/pollardslambda$ ./lambda <compressed public key(hex)> <key range(int)> <walkers(int)> <OPTIONAL DP(int)> <OPTIONAL Threads(int)> <OPTIONAL snaptime(int)>
    ```

    Replace `<compressed public key>` with the point \(G\) on the secp256k1 curve multiplied by your private key value, and `<key range>` with the size of the search interval for \(k\).

    Example usage:
    ```bash
    ~/pollardslambda$ ./lambda --pubkey 02145d2611c823a396ef6712ce0f712f09b9b4f3135e3e0aa3230fb9b6d08d1e16 --keyrange 135 --walkers 1000000 --dp 12 --t 8 --snaptime 15
    ```

## Commands

 The random walk begins using the public point of the compressed public key as the parameter H, the target private key range for initializing the initial probability space, and the optional distinguished points parameter, which will be calculated automatically if not defined:

```bash
~/pollardslambda$ ./lambda <compressed public key> <key range> <walkers> <dp bits> <threads> <snaptime>
```

--pubkey: The public key derived from the private key (discrete logarithm k) G = Q.

--keyrange: The range covering the regions of the elliptic curve where the discrete logarithm k resides.

--walkers: The walkers have the mission of traversing the elliptic curve point by point, they carry information such as the current point R and the coefficients a and b used to recover k in the event of a collision.

--dp: The Distinguished Points strategy is a memory-saving filter. Instead of storing every step of the walk (which would crash your RAM).

--t: Number of CPU threads/cores used to run the program.

--snaptime: Interval between each progress save, ```--snaptime 0``` don't save anything.

## External Libraries Used

"secp256k1.h" ```Lucas Leblanc```
"parallel_hashmap/phmap.h" ```Gregory Popovitch```

## Contributing

Contributions are welcome! Feel free to open issues or submit pull requests.

## Add a Star: <a href="https://github.com/lucaselblanc/pollardslambda/stargazers"><img src="https://img.shields.io/github/stars/lucaselblanc/pollardslambda?style=flat-square" alt="GitHub stars" style="vertical-align: bottom; width: 65px; height: auto;"></a>

## Donations: bc1pxqwuyfwvttjgttfmpt0gk0n7yzw3k7cyzzpc3rsc4lumr8ywythsj0rrhd

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

<p align="center">
  <a href="https://github.com/lucaselblanc">
    <img src="https://readme-typing-svg.demolab.com?font=Georgia&size=18&duration=2000&pause=100&multiline=true&width=500&height=80&lines=Lucas+Leblanc;Programmer+%7C+Student+%7C+Cyber+Security;+%7C+Android+%7C+Apps" alt="Typing SVG" />
  </a>
</p>

<a href="https://github.com/lucaselblanc">
    <img src="https://github-stats-alpha.vercel.app/api?username=lucaselblanc&cc=22272e&tc=37BCF6&ic=fff&bc=0000">
</a>

- 🔭 I’m currently working on [Pollard's Lambda Algorithm](https://github.com/lucaselblanc/pollardslambda/)

- 🚀 I’m looking to collaborate on: [Cyber-Security](https://play.google.com/store/apps/details?id=com.hashsuite.droid)

- 📝 I regularly read: [Monte Carlo methods for index computation (mod p)](https://www.ams.org/journals/mcom/1978-32-143/S0025-5718-1978-0491431-9/S0025-5718-1978-0491431-9.pdf)

- 📄 Know about my experiences: [https://www.linkedin.com/in/lucas-leblanc-215594208](https://www.linkedin.com/in/lucas-leblanc-215594208)

<br>

## GitHub Stats

<img src="https://streak-stats.demolab.com?user=lucaselblanc&theme=dracula"/>

<h3 align="left">Connect with me:</h3>
<p align="left">
<a href="https://www.linkedin.com/in/lucas-leblanc-215594208" target="blank"><img align="center" src="https://raw.githubusercontent.com/rahuldkjain/github-profile-readme-generator/master/src/images/icons/Social/linked-in-alt.svg" alt="lucas-leblanc-215594208" height="30" width="40" /></a>
<a href="https://www.youtube.com/@noclipstudiobr" target="blank"><img align="center" src="https://raw.githubusercontent.com/rahuldkjain/github-profile-readme-generator/master/src/images/icons/Social/youtube.svg" alt="@noclipstudiobr" height="30" width="40" /></a>
<a href="https://discord.gg/https://discord.gg/wXqcJDHht8" target="blank"><img align="center" src="https://raw.githubusercontent.com/rahuldkjain/github-profile-readme-generator/master/src/images/icons/Social/discord.svg" alt="https://discord.gg/wXqcJDHht8" height="30" width="40" /></a>
</p>

<h3 align="left">Languages and Tools:</h3>
<p align="left"> <a href="https://developer.android.com" target="_blank" rel="noreferrer"> <img src="https://raw.githubusercontent.com/devicons/devicon/master/icons/android/android-original-wordmark.svg" alt="android" width="40" height="40"/> </a> <a href="https://www.cprogramming.com/" target="_blank" rel="noreferrer"> <img src="https://raw.githubusercontent.com/devicons/devicon/master/icons/c/c-original.svg" alt="c" width="40" height="40"/> </a> <a href="https://www.w3schools.com/cpp/" target="_blank" rel="noreferrer"> <img src="https://raw.githubusercontent.com/devicons/devicon/master/icons/cplusplus/cplusplus-original.svg" alt="cplusplus" width="40" height="40"/> </a> <a href="https://www.w3schools.com/cs/" target="_blank" rel="noreferrer"> <img src="https://raw.githubusercontent.com/devicons/devicon/master/icons/csharp/csharp-original.svg" alt="csharp" width="40" height="40"/> </a> <a href="https://www.python.org" target="_blank" rel="noreferrer"> <img src="https://raw.githubusercontent.com/devicons/devicon/master/icons/python/python-original.svg" alt="python" width="40" height="40"/> </a> <a href="https://www.cprogramming.com/" target="_blank" rel="noreferrer"> <img src="https://raw.githubusercontent.com/devicons/devicon/master/icons/firebase/firebase-original.svg" alt="firebase" width="40" height="40"/> </a> <a href="https://www.linux.org/" target="_blank" rel="noreferrer"> <img src="https://raw.githubusercontent.com/devicons/devicon/master/icons/linux/linux-original.svg" alt="linux" width="40" height="40"/> </a> <a href="https://unity.com/" target="_blank" rel="noreferrer"> <img src="https://www.vectorlogo.zone/logos/unity3d/unity3d-icon.svg" alt="unity" width="40" height="40"/> </a> </p>