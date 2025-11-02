DotPlot(NKT, features = list(
    "T-like NKT" =
        c("CD3E","TRAC", "CD8A", "IL7R","KLF2"),
    "NK-like NKT" =
        c("NKG7", "GNLY", "PRF1","KLRD1","CCL5"),
    "Activated" = c("FOS", "JUN", "CD69"),
    "Exhausted" = c("PDCD1","LAG3","TIGIT")
)) + RotatedAxis()
