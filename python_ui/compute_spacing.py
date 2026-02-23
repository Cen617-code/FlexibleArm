import math

L0 = 55e-3
h0 = 350e-3

sqrt3 = math.sqrt(3.0)


def cosbeta1(h1, h2, h3, L0):
    delta12 = h2 - h1
    delta13 = h3 - h1

    p1p2 = (sqrt3 * L0 / 2.0, 3.0 * L0 / 2.0, delta12)
    xita1 = (0.0, sqrt3 * delta12 * L0, -1.5 * sqrt3 * L0 * L0)
    xita2 = (
        1.5 * delta13 * L0,
        sqrt3 * (delta12 * L0 - 0.5 * delta13 * L0),
        -1.5 * sqrt3 * L0 * L0,
    )

    crossVec = (
        xita1[1] * xita2[2] - xita1[2] * xita2[1],
        xita1[2] * xita2[0] - xita1[0] * xita2[2],
        xita1[0] * xita2[1] - xita1[1] * xita2[0],
    )
    dotValue = sum(a * b for a, b in zip(xita1, xita2))
    norm1 = math.sqrt(sum(a * a for a in xita1))
    norm2 = math.sqrt(sum(a * a for a in xita2))
    cosBeta = max(min(dotValue / (norm1 * norm2), 1.0), -1.0)
    beta = math.acos(cosBeta)
    triple_product = (
        p1p2[0] * crossVec[0]
        + p1p2[1] * crossVec[1]
        + p1p2[2] * crossVec[2]
    )
    if triple_product < 0:
        beta = -beta
    return beta


def niu(initialGuess, L0, L1, x13, y13):
    tolerance = 1e-20
    finiteStep = 1e-5
    maxIterations = 200

    x, y = initialGuess
    vectorP1P3 = (x13, y13)
    vectorP1P2 = (sqrt3 * L0 / 2.0, L0 / 2.0 + L1)

    def residual(point):
        px, py = point
        ocpP1 = (-sqrt3 * L0 / 2.0 - px, -L0 / 2.0 - py)
        termA = (ocpP1[0] + vectorP1P3[0], ocpP1[1] + vectorP1P3[1])
        termB = (ocpP1[0] + vectorP1P2[0], ocpP1[1] + vectorP1P2[1])
        normA = math.sqrt(termA[0] ** 2 + termA[1] ** 2)
        normB = math.sqrt(termB[0] ** 2 + termB[1] ** 2)
        normOc = math.sqrt(ocpP1[0] ** 2 + ocpP1[1] ** 2)
        residual1 = (
            0.5 * normA * normB
            + termA[0] * termB[0]
            + termA[1] * termB[1]
        )
        residual2 = (
            0.5 * normOc * normB
            + ocpP1[0] * termB[0]
            + ocpP1[1] * termB[1]
        )
        return (residual1, residual2)

    for _ in range(maxIterations):
        r1, r2 = residual((x, y))
        if r1 * r1 + r2 * r2 <= tolerance:
            break
        jac = [[0.0, 0.0], [0.0, 0.0]]
        for idx in range(2):
            if idx == 0:
                pt = (x + finiteStep, y)
            else:
                pt = (x, y + finiteStep)
            rp = residual(pt)
            jac[0][idx] = (rp[0] - r1) / finiteStep
            jac[1][idx] = (rp[1] - r2) / finiteStep
        det = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0]
        if abs(det) < 1e-10:
            # pseudo-inverse fallback could be added here if needed
            break
        inv_det = 1.0 / det
        dx = -(jac[1][1] * r1 - jac[0][1] * r2) * inv_det
        dy = -(-jac[1][0] * r1 + jac[0][0] * r2) * inv_det
        x += dx
        y += dy
    return x, y


def gammafast(h1, h2, h3, L0, L1):
    delta12 = h1 - h2
    sqrtA = math.sqrt(max(4.0 * delta12 * delta12 + 9.0 * L0 * L0, 0.0))
    sqrtB_arg = L0 * L0 * (
        4 * h1 * h1
        - 4 * h1 * h2
        - 4 * h1 * h3
        + 4 * h2 * h2
        - 4 * h2 * h3
        + 4 * h3 * h3
        + 9 * L0 * L0
    )
    sqrtB = math.sqrt(max(sqrtB_arg, 0.0))
    denominator = 4 * ((h1 - h2) ** 2 + 3 * L0 * L0) * L0
    term1 = 2 * h1 * h1 - 2 * h1 * h2 - 2 * h1 * h3 + 2 * h2 * h3 + 3 * L0 * L0
    numeratorX = term1 * L0 * L0 + sqrtA * sqrtB
    x13 = sqrt3 * numeratorX / denominator
    numeratorY = (
        2 * sqrtA * h1 * h1
        - 2 * sqrtA * h1 * h2
        - 2 * sqrtA * h1 * h3
        + 2 * sqrtA * h2 * h3
        + 3 * sqrtA * L0 * L0
        - 3 * sqrtB
    )
    denominatorY = 4 * ((h1 - h2) ** 2 + 3 * L0 * L0)
    y13 = numeratorY / denominatorY
    offsetX, offsetY = niu((0.0, 0.0), L0, L1, x13, y13)

    vectorOcP2 = (0.0, L1, 0.0)
    vectorOcpP2 = (-offsetX, L1 - offsetY, 0.0)
    cross_z = vectorOcP2[0] * vectorOcpP2[1] - vectorOcP2[1] * vectorOcpP2[0]
    dotValue = vectorOcP2[0] * vectorOcpP2[0] + vectorOcP2[1] * vectorOcpP2[1]
    normProduct = (
        math.sqrt(vectorOcP2[0] ** 2 + vectorOcP2[1] ** 2)
        * math.sqrt(vectorOcpP2[0] ** 2 + vectorOcpP2[1] ** 2)
    )
    cosGamma = max(min(dotValue / normProduct, 1.0), -1.0)
    gamma = math.acos(cosGamma)
    if cross_z < 0:
        gamma = -gamma
    return gamma, offsetX, offsetY


def jieA(h1, h2, h3, h0, L0):
    delta12 = h2 - h1
    halfBase = 1.5 * L0

    A1 = [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, h0],
        [0.0, 0.0, 0.0, 1.0],
    ]
    A2 = [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, h1 - h0],
        [0.0, 0.0, 0.0, 1.0],
    ]
    alpha = math.atan2(delta12, halfBase)
    cosA = math.cos(alpha)
    sinA = math.sin(alpha)
    A3 = [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, cosA, -sinA, -L0 / 2.0 + L0 * cosA / 2.0],
        [0.0, sinA, cosA, L0 * sinA / 2.0],
        [0.0, 0.0, 0.0, 1.0],
    ]

    L1 = math.sqrt(halfBase ** 2 + delta12 ** 2) - L0 / 2.0
    b = cosbeta1(h1, h2, h3, L0)

    triangleEdge = (sqrt3 * L0 / 2.0, L0 / 2.0 + L1, 0.0)
    edgeNorm = math.sqrt(sum(x * x for x in triangleEdge))
    unitP1 = triangleEdge[0] / edgeNorm
    unitP2 = triangleEdge[1] / edgeNorm

    denominator = math.sqrt(L0 * L0 + L1 * L1 + L0 * L1)
    Lq = (sqrt3 * L0 * L1 / 2.0) / denominator
    offsetX = -Lq * (1 - math.cos(b)) * ((L1 + L0 / 2.0) / denominator)
    offsetY = Lq * (1 - math.cos(b)) * ((sqrt3 * L0 / 2.0) / denominator)
    offsetZ = -Lq * math.sin(b)

    cosB = math.cos(b)
    sinB = math.sin(b)
    A4 = [
        [
            unitP1 * unitP1 * (1 - cosB) + cosB,
            unitP1 * unitP2 * (1 - cosB),
            unitP2 * sinB,
            offsetX,
        ],
        [
            unitP1 * unitP2 * (1 - cosB),
            unitP2 * unitP2 * (1 - cosB) + cosB,
            -unitP1 * sinB,
            offsetY,
        ],
        [-unitP2 * sinB, unitP1 * sinB, cosB, offsetZ],
        [0.0, 0.0, 0.0, 1.0],
    ]

    gamma, shiftX, shiftY = gammafast(h1, h2, h3, L0, L1)
    cosG = math.cos(gamma)
    sinG = math.sin(gamma)
    A5 = [
        [cosG, -sinG, 0.0, shiftX],
        [sinG, cosG, 0.0, shiftY],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]

    def matmul(a, b):
        result = [[0.0] * 4 for _ in range(4)]
        for i in range(4):
            for j in range(4):
                result[i][j] = sum(a[i][k] * b[k][j] for k in range(4))
        return result

    temp = matmul(A1, A2)
    temp = matmul(temp, A3)
    temp = matmul(temp, A4)
    temp = matmul(temp, A5)
    return temp


def rotation_z(angle):
    c = math.cos(angle)
    s = math.sin(angle)
    return [
        [c, -s, 0.0, 0.0],
        [s, c, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]


def matmul(a, b):
    result = [[0.0] * 4 for _ in range(4)]
    for i in range(4):
        for j in range(4):
            result[i][j] = sum(a[i][k] * b[k][j] for k in range(4))
    return result


def chain_transforms(stage_count, orientation_offset, h_init):
    positions = []
    transform = [[1.0 if i == j else 0.0 for j in range(4)] for i in range(4)]
    for stage in range(stage_count):
        orientation_angle = stage * orientation_offset
        local = jieA(h_init, h_init, h_init, h0, L0)
        rot = rotation_z(orientation_angle)
        stage_transform = matmul(rot, local)
        transform = matmul(transform, stage_transform)
        positions.append((transform[0][3], transform[1][3], transform[2][3]))
    return positions


if __name__ == "__main__":
    stage_count = 8
    orientation_offset = math.pi / 3
    positions = chain_transforms(stage_count, orientation_offset, h0)
    print("Stage center positions (m):")
    for idx, pos in enumerate(positions, 1):
        print(f"Stage {idx}: x={pos[0]:.6f}, y={pos[1]:.6f}, z={pos[2]:.6f}")
    print("\nAdjacent spacing (Euclidean distance, m):")
    prev = (0.0, 0.0, 0.0)
    for idx, pos in enumerate(positions, 1):
        dx = pos[0] - prev[0]
        dy = pos[1] - prev[1]
        dz = pos[2] - prev[2]
        dist = math.sqrt(dx * dx + dy * dy + dz * dz)
        print(f"Base to Stage {idx}: {dist:.6f}")
        prev = pos
    print("\nStage-to-stage spacing (m):")
    for idx in range(1, len(positions)):
        a = positions[idx - 1]
        b = positions[idx]
        dx = b[0] - a[0]
        dy = b[1] - a[1]
        dz = b[2] - a[2]
        dist = math.sqrt(dx * dx + dy * dy + dz * dz)
        print(f"Stage {idx} -> {idx+1}: {dist:.6f}")
