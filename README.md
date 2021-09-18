# Resistor network voltage divider optimization

These are some experiments with using discrete optimization techniques to a resitor network to approximate
a set of target potential levels. The design works around the assumption that a set of resistors will be
individually assigned to a pin of a microcontroller and the other ends all tied together. 
The microcontroller can selectively set each resistor pin to 5V (HIGH) or GND (LOW). This creates a 
voltage divider of two sets of parallel resistors. We can approximate up to 2<sup>N</sup> potentials 
with only N pins. However, the larger the set of desired potentials compared to N the less the overall
accuracy/convergence will be.
