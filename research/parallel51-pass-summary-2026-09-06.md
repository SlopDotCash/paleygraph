# Pass 51: collisions persist outside the newest primitive layer

The Paley conjecture and Proximity Prize remain unproved. This pass
rules out a shortcut in the new tower criterion and records an exact
constraint linking primitive labels at different subgroup sizes.

At p=2144280833, which lies in the order-256 quartic interval, the
order-128 subgroup has primitive triple mass 3, while the order-256
primitive polynomial has only simple roots. The full order-256 energy
retains the earlier excess: E_256=208128 and normalized excess D_256=48.
Thus total tower triple mass cannot be bounded by a constant times the
largest level's primitive triple mass. This does not contradict any
energy bound with a quadratic baseline or the full conjecture.

The [derivation](parallel51-collision-transport-2026-09-06.md) transports
each lower primitive label by the appropriate power to the endpoint.
It separates intrinsic collisions, mergers under those power maps,
coincidences across levels, and the distinguished label. An exact
compatibility equation shows these data cannot be assigned independently.
It does not yet provide a uniform estimate for them.

Retaining the opposite-pair baseline also sharpens the earlier lower
bound to

    E_N >= 3N^2-3N + 6N sum_s Y_s.

No uniform exponent improves. The unresolved task is still the upper
tower estimate from pass 50, now with an explicit reason not to replace
it by the newest primitive layer alone.

The [checker](../experiments/parallel51_collision_transport.py) verifies
71 field towers, 397 levels and 2,731,208 literal pairs. All checks pass,
including a separate primality certificate and complete triple witness.
The [results](../results/parallel51_collision_transport_2026_09_06.json)
and [audit](../results/parallel51_pass_audit_2026_09_06.json) preserve
the exact data and scope. The prime was already certified at order128;
its order256 endpoint is the new comparison, not a new classification.

Root completed ordinary proofs and checks. Previously usage-limited
agents were not restarted. No Lean process was polled or changed and
no Prove2Me submission was made. No independent-agent, Lean, external
review or novelty claim is made. The full goal remains active and unproved.
