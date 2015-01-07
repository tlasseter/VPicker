package com.Sts.SeismicAttributes;

import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.DataCube.*;
import static com.Sts.Utilities.DataCube.StsFileBlocks.*;
import com.Sts.Utilities.*;

import java.io.*;
import java.nio.*;
import java.nio.channels.*;

/**
 * Created by IntelliJ IDEA.
  * User: Tom Lasseter
  * Date: Jul 14, 2007
  * Time: 8:59:18 PM
  * To change this template use File | Settings | File Templates.
  */

/**
 * This class builds an StsSeismicVolume from input data and saves it in output files using mappedByteBuffers in the process.
  * Concrete subclasses implement allocateMemory and doProcessBlock methods to perform the work.
  * Input volumes and output volumes are the same size except in subclass StsPostStackSegyVolumeConstructor which can define
  * a smaller and/or decimated volume for output.
  */
 public abstract class StsVolumeConstructor
 {
     /** convenience copy of the model */
     StsModel model;
     /** name of the volume we are creating */
     String volumeName;
     /** volume we are creating */
     StsSeismicVolume outputVolume;

     /** row-ordered collection of file-mapped float blocks which is same as outputFloatBuffer */
     StsCubeFileBlocks rowFloatFileMapBlocks;
     /** row-ordered collection of file-mapped byte blocks which are written to on second pass from the float blocks */
     StsCubeFileBlocks rowByteFileMapBlocks;
     /** col-ordered collection of file-mapped byte blocks which are written to on second pass from the float blocks */
     StsCubeFileBlocks colByteFileMapBlocks;
     /** slice-ordered collection of file-mapped byte blocks which are written to on second pass from the float blocks */
     StsCubeFileBlocks sliceByteFileMapBlocks;
     /** buffer to which output float data is written */
     StsMappedFloatBuffer outputFloatBuffer;
     /** file to which col bytes are being written */
     RandomAccessFile colByteFile;
     /** file to which slice bytes are being written */
     RandomAccessFile sliceByteFile;

     boolean inputDataIsFloat = false;
     boolean outputDataIsFloat = true;
     StsProgressPanel panel = null;
     int nTracesDone;
     int nInputRows;
     int nInputCols;
     //    int inputRowInc;
     //    int inputColInc;
     //    int inputSliceInc;
     int nInputSlices;
     int nInputBlocks;
     int nInputRowsPerBlock;
     int nInputSamplesPerBlock;
     int nInputBlockTraces;
     int nInputBlockFirstRow;
     int nInputBlockLastRow;

     int nOutputRows;
     int nOutputCols;
     int nOutputSlices;

     int nOutputSamplesPerRow;
     int nOutputSamplesPerInputBlock;

     int nOutputBlocks;
     int nOutputBlockFirstRow;
     int nOutputBlockLastRow;
     int nOutputRowsPerBlock;
     int nOutputSamplesPerBlock;
     int nOutputBlockTraces;
     int nOutputRowPlaneSamples;
     int nOutputBlockColPlaneSamples;
     int nOutputBlockSlicePlaneSamples;
     //    int nOutputFloatBytesPerInputBlock;
     //    int nInputSamplesPerBlock;

     /** number of blocks in memory waiting to be written out for xlines and slices */
     int nBlocksInMemory;

     //    RandomAccessFile colByteFile = null;
     //    RandomAccessFile sliceByteFile = null;
     //   RandomAccessFile nativeRowFile = null;
     byte[] planeRowBytes = null;
     byte[] blockColBytes = null;
     byte[] blockSliceBytes = null;
     byte[] traceBytes = null;
     float[] traceFloats = null;

     long inputPosition;
     long outputPosition;

     float dataMin = StsParameters.largeFloat;
     float dataMax = -StsParameters.largeFloat;

     StsProgressDialog dialog;
     public boolean canceled = false;  // called to cancel the process

     // If appropriate, print out the block number as input processing runs
     protected boolean displayProcessingBlock = true;

     boolean debug = false;
     StsVolumeConstructorTimer timer;
     static public boolean runTimer = false;

     static final boolean usingOutputMappedBuffers = false;

     StsMemAllocVolumeProcess memoryAllocation;

     public abstract boolean allocateMemory();

     public abstract boolean processBlockInput(int nBlock, String mode);

     public StsSeismicVolume getVolume()
     {
         return outputVolume;
     }

     public abstract boolean isOutputDataFloat();

     public boolean initializeBlockInput(int nBlock)
     {
         nInputBlockFirstRow = nInputBlockLastRow + 1;
         nInputBlockLastRow += nInputRowsPerBlock;

         if (nInputBlockLastRow > nInputRows - 1)
         {
             nInputBlockLastRow = nInputRows - 1;
             nInputRowsPerBlock = nInputBlockLastRow - nInputBlockFirstRow + 1;
             nInputSamplesPerBlock = nInputRowsPerBlock * nInputCols * nInputSlices;
             nInputBlockTraces = nInputRowsPerBlock * nInputCols;
         }
         return true;
     }

     public void initializeBlockOutput(int nBlock)
     {
         nOutputBlockFirstRow = nOutputBlockLastRow + 1;
         nOutputBlockLastRow += nOutputRowsPerBlock;

         if (nOutputBlockLastRow > nOutputRows - 1)
         {
             nOutputBlockLastRow = nOutputRows - 1;
             nOutputRowsPerBlock = nOutputBlockLastRow - nOutputBlockFirstRow + 1;
             nOutputSamplesPerBlock = nOutputRowsPerBlock * nOutputRowPlaneSamples;
             nOutputBlockTraces = nOutputRowsPerBlock * nOutputCols;
             nOutputBlockColPlaneSamples = nOutputRowsPerBlock * nOutputSlices;
             nOutputBlockSlicePlaneSamples = nOutputRowsPerBlock * nOutputCols;
         }
     }

     public void createOutputVolume(String mode)
     {
		 final String mo = mode;
         StsProgressRunnable runnable = new StsProgressRunnable()
         {
             public void run()
             {
                 runCreateOutputVolume(mo);
             }

             public void cancel()
             {
                 panel.cancel();
             }
         };
         StsToolkit.runRunnable(runnable);
     }

	 public void createOutputVolume()
	 {
		 StsProgressRunnable runnable = new StsProgressRunnable()
		 {
			 public void run()
			 {
				 runCreateOutputVolume();
			 }

			 public void cancel()
			 {
				 panel.cancel();
			 }
		 };
		 StsToolkit.runRunnable(runnable);
	 }

	 protected boolean runCreateOutputVolume() { return runCreateOutputVolume(null); }

     protected boolean runCreateOutputVolume(String mode)
     {
         long startTime = 0, stopTime = 0;

         try
         {
             if (!allocateMemory()) return false;
             if (isCanceled()) return false;

             if (runTimer)
             {
                 StsVolumeConstructorTimer.clear();
                 StsVolumeConstructorTimer.overallTimer.start();
             }

             model.disableDisplay();

             outputVolume.clearHistogram();
             // In the first block loop (input), we compute a row float volume (if float data available) or a row byte volume (if not).
             // For a float volume, the doProcessBlock method should determine the min and max range for the volume
             // which can be subsequently used in scaling the row, col, and slice byte volumes produced in the second block loop (output).

             initializeVolumeInput();
             for (int nBlock = 0; nBlock < nInputBlocks; nBlock++)
             {
                 Main.logUsageTimer();
                 if (isCanceled()) return false;
                 initializeBlockInput(nBlock);
                 if (!processBlockInput(nBlock, mode)) return false;
                 if (panel != null)
                 {
                     if (displayProcessingBlock)
                         panel.appendLine("Processing input block " + (nBlock + 1) + " of " + nInputBlocks + ".");
                     panel.incrementCount();
                 }
                 finalizeBlockInput(nBlock);
             }
             finalizeVolumeInput();
             // we have computed a row float volume or a row byte volume; if a float row volume, we scale it and output row, col, and slice byte volumes;
             // if only a row byte volume, we output col and slice byte volumes
             initializeVolumeOutput();
             panel.appendLine("Output processing will be in " + nOutputBlocks + " blocks with " + nOutputRowsPerBlock + " rows per block.");
             //            panel.appendLine("    Writing " + memoryAllocation.nBlocksPerWrite + " input blocks per output block");
             panel.incrementInterval();
             panel.setIntervalCount(nOutputRows);
             for (int nBlock = 0; nBlock < nOutputBlocks; nBlock++)
             {
                 Main.logUsage();
                 if (isCanceled()) return false;
                 initializeBlockOutput(nBlock);
                 if (panel != null)
                     panel.appendLine("Writing block " + (nBlock + 1) + " of " + nOutputBlocks + " to disk.");
                 outputByteVolumes(nBlock);
                 finalizeBlockOutput(nBlock);
             }
             Main.logUsage();
             System.out.println("Module: " + Main.usageModule + " Message: " + Main.usageMessage);
             finalizeVolumeOutput();
             outputVolume.initialize();

             return true;
         }
         catch (Exception e)
         {
             StsException.outputWarningException(this, "initializeOutputVolume", e);
             String message = new String("Failed to process volume " + outputVolume.getName() + ".\n" + "Error: " + e.getMessage());
             new StsMessage(null, StsMessage.WARNING, message);
             return false;
         }
         finally
         {
             if (runTimer) timer.overallTimer.stopPrint("Total time to process seismic data: ");
             //            if(attributeVolume != null) attributeVolume.delete();
             if (dialog != null) dialog.dispose();
             StsVolumeConstructorTimer.printTimers("Volume: " + volumeName);
             model.enableDisplay();
         }
     }

     protected void initializeVolumeInput()
     {
         initializeGridRange();
         //        nInputRows = memoryAllocation.nInputRows;
         //        nInputCols = memoryAllocation.nInputCols;
         //        nInputSlices = memoryAllocation.nInputSlices;

         nInputBlockLastRow = -1;
         nInputRowsPerBlock = memoryAllocation.nInputRowsPerBlock;
         nInputBlocks = memoryAllocation.nInputBlocks;
         nInputBlockTraces = memoryAllocation.nInputBlockTraces;
         nInputSamplesPerBlock = memoryAllocation.nInputSamplesPerBlock;
         nOutputSamplesPerInputBlock = memoryAllocation.nOutputSamplesPerInputBlock;
         dataMin = StsParameters.largeFloat;
         dataMax = -StsParameters.largeFloat;
         traceFloats = new float[nInputSlices];

         panel.appendLine("Creating volume: " + outputVolume.getName());
         panel.appendLine("Volume size: nInputRows " + nInputRows + " nInputCols " + nInputCols + " nSamples " + nInputSlices);
         panel.appendLine("Input processing will be in " + nInputBlocks + " blocks with " + nInputRowsPerBlock + " rows per block.");
         panel.initializeIntervals(2, 100);
         panel.setIntervalCount(nInputBlocks);
     }

     protected void initializeGridRange()
     {
         nInputRows = memoryAllocation.nInputRows;
         nInputCols = memoryAllocation.nInputCols;
         nInputSlices = memoryAllocation.nInputSlices;
         nOutputRows = nInputRows;
         nOutputCols = nInputCols;
         nOutputSlices = nInputSlices;
         //        inputRowInc = 1;
         //        inputColInc = 1;
         //        inputSliceInc = 1;
     }

     protected void initializeVolumeOutput()
     {
         nOutputBlockLastRow = -1;
         outputVolume.checkScaling();
         nOutputRowsPerBlock = memoryAllocation.nOutputRowsPerBlock;
         nOutputBlocks = memoryAllocation.nOutputBlocks;
         nOutputRowPlaneSamples = nOutputCols * nOutputSlices;
         nOutputSamplesPerBlock = nOutputRowsPerBlock * nOutputRowPlaneSamples;
         nOutputBlockColPlaneSamples = nOutputRowsPerBlock * nOutputSlices;
         nOutputBlockSlicePlaneSamples = nOutputRowsPerBlock * nOutputCols;
         rowFloatFileMapBlocks = outputVolume.fileMapRowFloatBlocks;
         StsCubeFileBlocks[] filesMapBlocks = outputVolume.filesMapBlocks;
         rowByteFileMapBlocks = filesMapBlocks[1];
         if (usingOutputMappedBuffers)
         {
             colByteFileMapBlocks = filesMapBlocks[0];
             sliceByteFileMapBlocks = filesMapBlocks[2];
         }
         else
         {
             colByteFile = filesMapBlocks[0].getFile();
             sliceByteFile = filesMapBlocks[2].getFile();
         }
         traceFloats = new float[nOutputSlices];
         traceBytes = new byte[nOutputSlices];
         planeRowBytes = new byte[nOutputRowPlaneSamples];
         blockColBytes = new byte[nOutputSamplesPerBlock];
         blockSliceBytes = new byte[nOutputSamplesPerBlock];
         outputVolume.initializeScaling();

     }

     public void adjustOutputDataRange(final float[] floatData)
     {
         for (int n = 0; n < floatData.length; n++)
         {
             float value = floatData[n];
             if (value < dataMin)
                 dataMin = value;
             else if (value > dataMax)
                 dataMax = value;
         }
     }

     private boolean outputByteVolumes(int nBlock)
     {
         int row, col, slice, blockRow, nBlockColOffset, nBlockSliceOffset;
         int n = 0;
         try
         {
             if (runTimer) StsVolumeConstructorTimer.getFloatRowTimer.start();

             for (row = nOutputBlockFirstRow, blockRow = 0; row <= nOutputBlockLastRow; row++, blockRow++)
             {
                 ByteBuffer outputRowByteBuffer;

                 if (rowFloatFileMapBlocks != null)
                 {
                     FloatBuffer outputRowFloatBuffer = rowFloatFileMapBlocks.getByteBufferPlane(row, FileChannel.MapMode.READ_ONLY).asFloatBuffer();
                     outputRowByteBuffer = rowByteFileMapBlocks.getByteBufferPlane(row, FileChannel.MapMode.READ_WRITE);
                     outputVolume.scaleFloatsToBytes(outputRowFloatBuffer, outputRowByteBuffer, nOutputSamplesPerRow);
                     outputRowByteBuffer.rewind();
                 }
                 else
                     outputRowByteBuffer = rowByteFileMapBlocks.getByteBufferPlane(row, FileChannel.MapMode.READ_ONLY);

                 for (col = 0; col < nOutputCols; col++)
                 {
                     outputRowByteBuffer.get(traceBytes);
                     System.arraycopy(traceBytes, 0, planeRowBytes, col * nOutputSlices, nOutputSlices);
                     System.arraycopy(traceBytes, 0, blockColBytes, col * nOutputBlockColPlaneSamples + blockRow * nOutputSlices, nOutputSlices);

                     n = blockRow * nOutputCols + col;
                     for (int k = 0; k < nOutputSlices; k++, n += nOutputBlockSlicePlaneSamples)
                     {
                         blockSliceBytes[n] = traceBytes[k];
                         outputVolume.accumulateHistogram(traceBytes[k]);
                     }
                 }
                 rowFloatFileMapBlocks.unlockPlane(row);
                 rowByteFileMapBlocks.unlockPlane(row);
                 panel.incrementCount();
             }
             if (runTimer)
                 StsVolumeConstructorTimer.getFloatRowTimer.stopAccumulateIncrementCountPrintInterval("Block " + nBlock + " Mapped buffer reading row floats, scaling & writing in memory byte blocks:");

             if (runTimer) StsVolumeConstructorTimer.getByteRowTimer.start();

             if (!usingOutputMappedBuffers)
             {
                 int nColPlaneSamples = nOutputRows * nOutputSlices;
                 for (col = 0, nBlockColOffset = 0; col < nOutputCols; col++, nBlockColOffset += nOutputBlockColPlaneSamples)
                 {
                     long position = (long) col * (long) nColPlaneSamples + (long) nOutputBlockFirstRow * (long) nOutputSlices;
                     colByteFile.seek(position);
                     colByteFile.write(blockColBytes, nBlockColOffset, nOutputBlockColPlaneSamples);
                 }

                 int nSamplePlaneSamples = nOutputRows * nOutputCols;
                 for (slice = 0, nBlockSliceOffset = 0; slice < nOutputSlices; slice++, nBlockSliceOffset += nOutputBlockSlicePlaneSamples)
                 {
                     long position = (long) slice * (long) nSamplePlaneSamples + (long) nOutputBlockFirstRow * (long) nOutputCols;
                     sliceByteFile.seek(position);
                     sliceByteFile.write(blockSliceBytes, nBlockSliceOffset, nOutputBlockSlicePlaneSamples);
                 }
                 if (runTimer)
                     StsVolumeConstructorTimer.getByteRowTimer.stopAccumulateIncrementCountPrintInterval("Block " + nBlock + " RAF output xline & slice blocks:");
             }
             else
             {
                 for (col = 0, nBlockColOffset = 0; col < nOutputCols; col++, nBlockColOffset += nOutputBlockColPlaneSamples)
                 {
                     colByteFileMapBlocks.putColCubePlaneValues(nOutputBlockFirstRow, col, nOutputBlockColPlaneSamples, nBlockColOffset, blockColBytes);
                     colByteFileMapBlocks.unlockPlane(col);
                 }
                 for (slice = 0, nBlockSliceOffset = 0; slice < nOutputSlices; slice++, nBlockSliceOffset += nOutputBlockSlicePlaneSamples)
                 {
                     sliceByteFileMapBlocks.putSliceCubePlaneValues(nOutputBlockFirstRow, slice, nOutputBlockSlicePlaneSamples, nBlockSliceOffset, blockSliceBytes);
                     sliceByteFileMapBlocks.unlockPlane(slice);
                 }
                 if (runTimer)
                     StsVolumeConstructorTimer.getByteRowTimer.stopAccumulateIncrementCountPrintInterval("Block " + nBlock + " Mapped buffer output xline & slice blocks:");
             }
             return true;
         }
         catch (Exception e)
         {
             StsException.outputWarningException(this, "outputByteVolumes", e);
             if (runTimer)
                 StsVolumeConstructorTimer.getByteRowTimer.stopAccumulateIncrementCountPrintInterval("Block " + nBlock + " Exception termination.");
             return false;
         }
     }

     private boolean outputByteVolumesOld(int nBlock)
     {
         int row, col, slice, blockRow, nBlockColOffset, nBlockSliceOffset;
         int n = 0;
         try
         {
             if (runTimer) StsVolumeConstructorTimer.getFloatRowTimer.start();

             if (rowFloatFileMapBlocks != null)
             {
                 for (row = nOutputBlockFirstRow, blockRow = 0; row <= nOutputBlockLastRow; row++, blockRow++)
                 {
                     FloatBuffer floatBuffer = rowFloatFileMapBlocks.getByteBufferPlane(row, FileChannel.MapMode.READ_ONLY).asFloatBuffer();

                     for (col = 0; col < nOutputCols; col++)
                     {
                         try
                         {
                             floatBuffer.get(traceFloats);
                         }
                         catch (Exception e)
                         {
                             int correctPosition = col * nOutputSlices;
                             int bufferPosition = floatBuffer.position();
                             int correctCapacity = nOutputCols * nOutputSlices;
                             int bufferCapacity = floatBuffer.capacity();
                             StsException.outputWarningException(this, "outputByteVolumes", " get(traceFloats) row: " + row + " col: " + col +
                                 " correct position: " + correctPosition + " bufferPosition: " + bufferPosition +
                                 " correct capacity: " + correctCapacity + " bufferCapacity: " + bufferCapacity, e);
                             return false;
                         }
                         outputVolume.scaleFloatsToBytes(traceFloats, traceBytes);
                         System.arraycopy(traceBytes, 0, planeRowBytes, col * nOutputSlices, nOutputSlices);
                         System.arraycopy(traceBytes, 0, blockColBytes, col * nOutputBlockColPlaneSamples + blockRow * nOutputSlices, nOutputSlices);

                         n = blockRow * nOutputCols + col;
                         for (int k = 0; k < nOutputSlices; k++, n += nOutputBlockSlicePlaneSamples)
                         {
                             blockSliceBytes[n] = traceBytes[k];
                             outputVolume.accumulateHistogram(traceBytes[k]);
                         }
                     }
                     rowByteFileMapBlocks.putRowCubePlaneValues(row, nOutputRowPlaneSamples, 0, planeRowBytes);
                     rowByteFileMapBlocks.unlockPlane(row);
                     panel.incrementCount();
                 }
             }
             else
             {
                 for (row = nOutputBlockFirstRow, blockRow = 0; row <= nOutputBlockLastRow; row++, blockRow++)
                 {
                     ByteBuffer outputRowByteBuffer = rowByteFileMapBlocks.getByteBufferPlane(row, FileChannel.MapMode.READ_ONLY);
                     for (col = 0; col < nOutputCols; col++)
                     {

                         try
                         {
                             outputRowByteBuffer.get(traceBytes);
                         }
                         catch (Exception e)
                         {
                             int correctPosition = col * nOutputSlices;
                             int bufferPosition = outputRowByteBuffer.position();
                             int correctCapacity = nOutputCols * nOutputSlices * 4;
                             int bufferCapacity = outputRowByteBuffer.capacity();
                             StsException.outputWarningException(this, "outputByteVolumes", "get(traceBytes) row: " + row + " col: " + col +
                                 " correct position: " + correctPosition + " bufferPosition: " + bufferPosition +
                                 " correct capacity: " + correctCapacity + " bufferCapacity: " + bufferCapacity, e);
                             return false;
                         }
                         System.arraycopy(traceBytes, 0, blockColBytes, col * nOutputBlockColPlaneSamples + blockRow * nOutputSlices, nOutputSlices);

                         n = blockRow * nOutputCols + col;
                         for (int k = 0; k < nOutputSlices; k++, n += nOutputBlockSlicePlaneSamples)
                         {
                             blockSliceBytes[n] = traceBytes[k];
                         }
                     }
                     panel.incrementCount();
                 }
             }
             if (runTimer)
                 StsVolumeConstructorTimer.getFloatRowTimer.stopAccumulateIncrementCountPrintInterval("Block " + nBlock + " Mapped buffer reading row floats, scaling & writing in memory byte blocks:");

             if (runTimer) StsVolumeConstructorTimer.getByteRowTimer.start();

             if (!usingOutputMappedBuffers)
             {
                 int nColPlaneSamples = nOutputRows * nOutputSlices;
                 for (col = 0, nBlockColOffset = 0; col < nOutputCols; col++, nBlockColOffset += nOutputBlockColPlaneSamples)
                 {
                     long position = (long) col * (long) nColPlaneSamples + (long) nOutputBlockFirstRow * (long) nOutputSlices;
                     colByteFile.seek(position);
                     colByteFile.write(blockColBytes, nBlockColOffset, nOutputBlockColPlaneSamples);
                 }

                 int nSamplePlaneSamples = nOutputRows * nOutputCols;
                 for (slice = 0, nBlockSliceOffset = 0; slice < nOutputSlices; slice++, nBlockSliceOffset += nOutputBlockSlicePlaneSamples)
                 {
                     long position = (long) slice * (long) nSamplePlaneSamples + (long) nOutputBlockFirstRow * (long) nOutputCols;
                     sliceByteFile.seek(position);
                     sliceByteFile.write(blockSliceBytes, nBlockSliceOffset, nOutputBlockSlicePlaneSamples);
                 }
                 if (runTimer) StsVolumeConstructorTimer.getByteRowTimer.stopAccumulateIncrementCountPrintInterval("Block " + nBlock + " RAF output xline & slice blocks:");
             }
             else
             {
                 for (col = 0, nBlockColOffset = 0; col < nOutputCols; col++, nBlockColOffset += nOutputBlockColPlaneSamples)
                 {
                     colByteFileMapBlocks.putColCubePlaneValues(nOutputBlockFirstRow, col, nOutputBlockColPlaneSamples, nBlockColOffset, blockColBytes);
                     colByteFileMapBlocks.unlockPlane(col);
                 }
                 for (slice = 0, nBlockSliceOffset = 0; slice < nOutputSlices; slice++, nBlockSliceOffset += nOutputBlockSlicePlaneSamples)
                 {
                     sliceByteFileMapBlocks.putSliceCubePlaneValues(nOutputBlockFirstRow, slice, nOutputBlockSlicePlaneSamples, nBlockSliceOffset, blockSliceBytes);
                     sliceByteFileMapBlocks.unlockPlane(slice);
                 }
                 if (runTimer) StsVolumeConstructorTimer.getByteRowTimer.stopAccumulateIncrementCountPrintInterval("Block " + nBlock + " Mapped buffer output xline & slice blocks:");
             }
             return true;
         }
         catch (Exception e)
         {
             StsException.outputWarningException(this, "outputByteVolumes", e);
             if (runTimer) StsVolumeConstructorTimer.getByteRowTimer.stopAccumulateIncrementCountPrintInterval("Block " + nBlock + " Exception termination.");
             return false;
         }
     }

     protected boolean isCanceled()
     {
         if (panel == null || !panel.isCanceled()) return false;
         blocksMemoryManager.clearAllBlocks();
         return true;
     }

     public void finalizeBlockInput(int nBlock)
     {
         if (runTimer) StsVolumeConstructorTimer.getClearInputBlockTimer.start();
         //        if (rowFileMapFloatBlocks != null) rowFileMapFloatBlocks.clearBlock(nBlock);
         //        if (rowFileMapFloatBlocks != null) rowFileMapByteBlocks.clearBlock(nBlock);
         blocksMemoryManager.clearUnlockedBlocks();
         //       if(outputRowFloatBuffer != null) cleanClearBuffer(outputRowFloatBuffer);
         //       if(outputRowByteBuffer != null) cleanClearBuffer(outputRowByteBuffer);
         if (runTimer) StsVolumeConstructorTimer.getClearInputBlockTimer.stopAccumulateIncrementCountPrintInterval("clearUnlockBlocks");
     }

     public void finalizeBlockOutput(int nBlock)
     {
         if (runTimer) StsVolumeConstructorTimer.getClearOutputBlockTimer.start();
         if (rowFloatFileMapBlocks != null) rowFloatFileMapBlocks.unlockBlock(nBlock);
         if (rowByteFileMapBlocks != null) rowByteFileMapBlocks.unlockClearAllBlocks();

         if (usingOutputMappedBuffers)
         {
             if (colByteFileMapBlocks != null) colByteFileMapBlocks.unlockAllBlocks();
             if (sliceByteFileMapBlocks != null) sliceByteFileMapBlocks.unlockAllBlocks();
         }
         if (runTimer) StsVolumeConstructorTimer.getClearOutputBlockTimer.stopAccumulateIncrementCountPrintInterval("clearUnlockAllBlocks");
     }

     protected void cleanClearBuffer(MappedByteBuffer buffer)
     {
         if (buffer == null) return;
         //       buffer.force();
         buffer.clear();
         StsToolkit.clean(buffer);
     }

     public void finalizeVolumeInput()
     {
     }

     public void finalizeVolumeOutput()
     {
         if (runTimer) StsVolumeConstructorTimer.getClearOutputBlockTimer.start();
         finalizeVolumeOutputBuffers();
         if (runTimer) StsVolumeConstructorTimer.getClearOutputBlockTimer.stopAccumulateIncrementCountPrintInterval("final output block clear");
         outputVolume.calculateHistogram();
         outputVolume.writeHeaderFile();
         if (!outputVolume.addToProject(true)) return;
         outputVolume.addToModel();
         model.getProject().runCompleteLoading();
     }

     public void finalizeVolumeOutputBuffers()
     {
         blocksMemoryManager.unlockAndClearBlocks(rowByteFileMapBlocks);
         if (usingOutputMappedBuffers)
         {
             blocksMemoryManager.unlockAndClearBlocks(colByteFileMapBlocks);
             blocksMemoryManager.unlockAndClearBlocks(sliceByteFileMapBlocks);
         }
         else
         {
             /*
             try
             {
                 colByteFileMapBlocks.closeFile();
                 sliceByteFileMapBlocks.closeFile();
             }
             catch(Exception e)
             {
                 StsException.systemError(this, "finalizeVolumeOutputBuffers", "Failed try to close files for col and slice byte files.");
             }
             */
         }
         blocksMemoryManager.clearUnlockedBlocks();
         planeRowBytes = null;
         blockColBytes = null;
         blockSliceBytes = null;
     }
 }
