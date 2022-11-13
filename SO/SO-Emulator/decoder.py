
args = {
    0: "A",
    1: "D",
    2: "X",
    3: "Y",
    4: "[X]",
    5: "[Y]",
    6: "[X + D]",
    7: "[Y + D]",
}

g_0_codes = {
    0x0000: "MOV",
    0x0002: "OR",
    0x0004: "ADD",
    0x0005: "SUB",
    0x0006: "ADC",
    0x0007: "SBB",
    0x0008: "XCHG"
}

def g_0(x):
    arg1 = (x & 0x0700) >> 8
    arg2 = (x & 0x3800) >> 11
    opcode = x & 0x00ff
    op = g_0_codes[opcode]
    
    return f"{op} {args[arg1]} {args[arg2]}"

g_1_codes = {
    0: "MOVI",
    3: "XORI",
    4: "ADDI",
    5: "CMPI",
    6: "RCR",
}

def g_1(x):
    arg1 = (x & 0x0700) >> 8
    imm8 = (x & 0x00ff)
    opcode = (x & 0x3800) >> 11
    op = g_1_codes[opcode]

    return f"{op} {args[arg1]} {imm8}"

def g_2(x):
    if x == 0x8000:
        return "CLC"
    elif x == 0x8100:
        return "STC"
    return ""

g_3_codes = {
    0: "JMP",
    2: "JNC",
    3: "JC",
    4: "JNZ",
    5: "JZ",
}

def g_3(x):
    imm8 = x & 0x00ff
    opcode = (x & 0x3f00) >> 8
    return f"{g_3_codes[opcode]} {imm8}"

def decode(ins):
    res = ""
    if ins & 0xc000 == 0x0000:
        res = g_0(ins)
    elif ins & 0xc000 == 0x4000:
        res = g_1(ins)
    elif ins & 0xc000 == 0x8000:
        res = g_2(ins)
    elif ins & 0xc000 == 0xc000:
        res = g_3(ins)
    print(res)


print(decode(0x0002 + 0x100 * 5 + 0x0800 * 3))
