package org.benchmark.filter4;

import javax.vecmath.Matrix4d;

public record CrystalNeighborHit(
        CrystalSearchTarget target,
        CrystalSourceObject source,
        CrystalInstanceKey key,
        String symmetryOperationString,
        Matrix4d transform,
        double minDistanceLigand,
        double minDistancePocket,
        int ligandContactCount,
        int pocketContactCount) {
}
