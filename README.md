Here is code for analyses of ovary ATAC data.

Following pipeline was used:
1. Combine samples and calculate gene scores using snapATAC2-based pipeline: https://github.com/cellgeni/nf-atac/tree/snapatac2
2. Combine with single cell reference
3. Use scanvi for integration and label transfer with dataset_donor as batch
4. For some reason scanvi label transfer misannotated many reference cells (for examples it annotated all Stroma_InnCor as Immune_Mac_resident_SPP1+) so I used bbknn on scanvi embedding to transfer labels majority voting of 50 nearest reference neighbors. Integration generally looks ok, but Granulosa cells, specifically AMH-expressed ones are almost complitely absent in ATAC.
