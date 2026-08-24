package org.benchmark.filter4;

public record CrystalInstanceKey(
        String sourceObjectId, int symmetryOperationId, int h, int k, int l) {
    @Override
    public String toString() {
        return sourceObjectId + "|op_" + symmetryOperationId + "|" +
                String.format("%+d|%+d|%+d", h, k, l);
    }
}
