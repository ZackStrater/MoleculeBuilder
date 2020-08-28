bond_encoder = {
    ".": 0,  # non-bond
    "=": 2,  # double bond
    "#": 3,  # triple bond
    # aromatic bonds (4) denoted by lowercase element symbol
    "$": 5,  # ionic bond
    "/": 6,  # vinyl up
    "\\": 7,  # vinyl down
}


def encode_bond(bonding_info):
    detected = False
    code = None
    for key in bond_encoder:
        if key in bonding_info:
            code = (bond_encoder[key])
            detected = True
    if not detected:
        return 1
    else:
        return code

print(encode_bond("."))