# Num2Bits

```
template Num2Bits (nBits) {
    signal input a;
    signal output b[nBits];

    var sum = 0;
    for (var i = 0; i< nBits; i++){
        b[i] <--  (a >> i) & 1;
        b[i] * (b[i] -1) === 0;
        sum = sum + 2**i*b[i];
    }

    a === sum;
}

component main { public [a]}  = Num2Bits(5);
/* INPUT = {
    "a": "5"
} */
```

# IsZero
```
template IsZero () {
    signal input in;
    signal output out;
    signal inv;
    inv <-- in != 0 ? 0 : 1;
    out <== inv;
    out*(out -1) === 0;
}

component main { public [in]}  = IsZero();

/* INPUT = {
    "in": "0"
} */
```

# IsEqual

```
template IsEqual () {
    signal input in[2];
    signal output out;
    component isz = IsZero();
    isz.in <== in[0] - in[1];
    out <== isz.out;
    out*(out -1) === 0;
}

component main { public [in]}  = IsEqual();

/* INPUT = {
    "in": ["2","0"]
} */
```

# Selector
```
template Selector(nChoices) {
    signal input in[nChoices];
    signal input index;

    signal output out;
    signal res;
    res <-- index < nChoices ? in[index] : 0;
    out <== res;
}

component main { public [in]}  = Selector(3);

/* INPUT = {
    "in": ["1","2","3"],
    "index": "2"
} */
```

# IsNegative
Mental model of (**negative**) integer denotation:
```
p = 21888242871839275222246405745257275088548364400416034343698204186575808495617.
field element range: 0-------p/2---------p
positive: 0, 1, ..., p/2-1
negative: p-1 => -1, p-2=>-2, p-3==> -3, ..., p-(p/2-1) => -(p/2 -1)
```
> Notice: In circom, assert(-1 == 21888242871839275222246405745257275088548364400416034343698204186575808495616) is Ok. However you can't just input `-1` in the input field of **/*INPUT** when you do testing, due to it's actually 4294967295 if you log the corresponding input.
```
template IsNegative(){
    signal input in;
    signal output out;

    component n2b = Num2Bits(254);
    component compc = CompConstant(10944121435919637611123202872628637544274182200208017171849102093287904247808);
    n2b.in <== in;
    compc.in <== n2b.out;
    out <== compc.out;


}

component main { public [in]}  = IsNegative();

/* INPUT = {
    "in": "-1"
} */
```

Q: Understanding check: Why can’t we just use LessThan or one of the comparator circuits from the previous exercise?

A: The defination of [Relational operator](https://docs.circom.io/circom-language/basic-operators/#boolean-operators) **<** constrains val(z) to be not greater than p/2. Thus, if we compare `x` (less than p/2) with `y` (larger than p/2) when use Less than, what Circom does is comparing `x` with `p - y`, then maybe we will get the wrong answer.

# LessThan

## Extension 1
### LessThan
```
template LessThan(k) {
    assert(k<=252);
    signal input in[2];
    signal output out;
    signal temp;
    temp <-- in[0] < in[1] ?  1 : 0;
    out <== temp;
}

component main { public [in]}  = LessThan(252);

/* INPUT = {
    "in": ["7237005577332262213973186563042994240829374041602535252466099000494570602493","7237005577332262213973186563042994240829374041602535252466099000494570602492"]
} */
```

## Extension 2
### LessEqThan
```
template LessEqThan() {
    signal input in[2];
    signal output out;

    component n2b = Num2Bits(253);

    n2b.in <== in[1]+ (1<<252) - in[0];

    out <==  n2b.out[252];
}

component main { public [in]}  = LessEqThan();

/* INPUT = {
    "in": ["7237005577332262213973186563042994240829374041602535252466099000494570602491","7237005577332262213973186563042994240829374041602535252466099000494570602492"]
} */
```

### GreaterThan
```
template GreaterThan() {
    signal input in[2];
    signal output out;

    component n2b = Num2Bits(253);

    n2b.in <== in[1]+ (1<<252) - in[0];

    out <== 1- n2b.out[252];
}

component main { public [in]}  = GreaterThan();

/* INPUT = {
    "in": ["7237005577332262213973186563042994240829374041602535252466099000494570602493","7237005577332262213973186563042994240829374041602535252466099000494570602492"]
} */
```

### GreaterEqThan
```
template GreaterEqThan() {
    signal input in[2];
    signal output out;

    component n2b = Num2Bits(253);

    n2b.in <== in[0]+ (1<<252) - in[1];

    out <== n2b.out[252];
}

component main { public [in]}  = GreaterEqThan();

/* INPUT = {
    "in": ["7237005577332262213973186563042994240829374041602535252466099000494570602492","7237005577332262213973186563042994240829374041602535252466099000494570602492"]
} */
```


# IntegerDivide

```
template IntegerDivide(divisor_bits) {
    assert(divisor_bits <= 126);

    signal input dividend; // -8
    signal input divisor; // 5
    signal output remainder; // 2
    signal output quotient; // -2

    component is_neg_comp = IsNegative();
    is_neg_comp.in <== dividend;
    var is_neg = is_neg_comp.out;

    signal output sign;
    sign <== 1 + is_neg * -2; // 1 or -1

    var abs_dividend = dividend * sign; // 8
    var raw_remainder = abs_dividend % divisor;
    var neg_remainder =  divisor - raw_remainder;

    remainder <-- is_neg == 1 && raw_remainder != 0 ? neg_remainder: raw_remainder;
    quotient <-- (dividend - remainder) / divisor; // (-8 - 2) / 5 = -2.
    dividend === divisor * quotient + remainder; // -8 = 5 * -2 + 2.

    // check that 0 <= remainder < divisor
    component remainderUpper = LessThan(divisor_bits);
    remainderUpper.in[0] <== remainder;
    remainderUpper.in[1] <== divisor;
    remainderUpper.out === 1;
}

component main { public [dividend, divisor]}  = IntegerDivide(8);
//-8: 21888242871839275222246405745257275088548364400416034343698204186575808495609
//-2: 21888242871839275222246405745257275088548364400416034343698204186575808495615
/* INPUT = {
    "dividend": "21888242871839275222246405745257275088548364400416034343698204186575808495609",
    "divisor": "5"
} */
```