package com.Sts.DBTypes;

import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.io.*;

/**
 * <p>Title: </p>
 * <p/>
 * <p>Description: </p>
 * <p/>
 * <p>Copyright: Copyright (c) 2007</p>
 * <p/>
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */

public class StsSeismicExportPanel3d extends StsSeismicExportPanel
{
    public StsSeismicExportPanel3d()
    {
        super();
    }

    public StsSeismicExportPanel3d(StsModel model, StsSeismicBoundingBox volume, String title, boolean isAVelocityExport)
    {
        super(model, volume, title, isAVelocityExport);
    }

    protected StsGroupBox getInputGeometryGroupBox()
    {
        return new StsSeismicVolume3dGroupBox(volume, "Selected Volume");
    }

    protected StsGroupBox getOutputGeometryGroupBox()
    {
        return new StsSeismicVolume3dGroupBox(exportVolume, "Export Volume", volume);
    }

    protected void doExportSeismic(StsProgressPanel progressPanel)
    {
        if(setupExport())
        {
            int progressRow = 0;
            StsSEGYFormat segyFmt = StsSEGYFormat.constructor(model, StsSEGYFormat.POSTSTACK);
            StsProject project = model.getProject();
            String filename = project.getRootDirString() + exportName + ".sgy";

            byte currentType = velocityType;
            // type is not saved in segy, so we need to find a place for it in the segy file headers
            volume.setType(StsParameters.SAMPLE_TYPE_VEL_RMS);
            initializeRowColExportRanges();
            try
            {
                progressPanel.setMaximum(croppedBoundingBox.nRows);
                progressPanel.appendLine("Outputting " + filename);

                setupFile(segyFmt, filename);
                progressPanel.appendLine("    Outputting text header.");
                writeTextHdr(segyFmt);
                progressPanel.appendLine("    Outputting binary header.");
                writeBinaryHeader();
                progressPanel.appendLine("    Outputting trace data." + filename);

                int rowMin = croppedBoundingBox.rowMin;
                int rowMax = croppedBoundingBox.rowMax;
                int rowInc = croppedBoundingBox.rowInc;
                float rowNum = croppedBoundingBox.rowNumMin;
                float rowNumInc = croppedBoundingBox.rowNumInc;
                int row = rowMin;
                int rowProgress = 0;
                while (!canceled)
                {
                    if(!exportRow(segyFmt, row, rowNum))
                    {
                        progressPanel.appendErrorLine("Failed to read row " + row + " inline " + rowNum + ".");
                        break;
                    }
                    ;
                    progressPanel.setValue(++rowProgress);
                    String msg = "Completed row " + row + " of " + rowMin + "-" + rowMax;
                    progressPanel.setDescription(msg);
                    if(row == rowMax)
                        break;
                    row += rowInc;
                    rowNum += rowNumInc;
                }
                progressPanel.appendLine("Completed output of sample data.");
            }
            catch(Exception e)
            {
                StsException.outputException("StsSeismic.export() failed.", e, StsException.WARNING);
                return;
            }
            finally
            {
                try
                {
                    dos.flush();
                    dos.close();
                    dos = null;
                }
                catch(Exception e)
                {
                    StsException.outputException("StsSeismicVolume.export() failed.", e, StsException.WARNING);
                    return;
                }
                progressPanel.setValue(nRows);
                out = StsFile.constructor(filename);
                long fileLength = out.length();
                progressPanel.appendLine("Completed: " + filename + " output.");
                progressPanel.appendLine("File Size: " + fileLength);

                long correctLength = textHeaderSize + binaryHeaderSize +
                        nRows * nCols * (traceHeaderSize + nCroppedSlices * sampleSize);
                System.out.println("fileLength: " + fileLength + " correctLength: " + correctLength);
            }
            // type is not saved in segy, so we need to find a place for it in the segy file headers
            volume.setType(currentType);
        }
    }

    private void initializeRowColExportRanges()
    {

    }

    protected boolean checkExportSeismic()
    {
        StsProject project = model.getProject();
        String filename = project.getRootDirString() + exportName + ".sgy";
        File file = new File(filename);
        if(file.exists())
        {
            boolean overWrite = StsYesNoDialog.questionValue(model.win3d, "File " + filename + " already exists. Do you wish to overwrite it?");
            if(!overWrite)
                return false;
            file.delete();
        }

        return true;
    }

    protected boolean setupExport()
    {
        computeConversionFlag();
        setGridParameters();
        constructTextHeader();
        constructBinaryHdr();

        isInputFloats = volume.isDataFloat;
        if(isInputFloats)
        {
            if(!volume.setupReadRowFloatBlocks())
            {
                if(!volume.setupRowByteBlocks())
                {
                    new StsMessage(model.win3d, StsMessage.WARNING,
                            "Failed to find float or byte data volume to read: cannot process.");
                    return false;
                }
                else
                {
                    new StsMessage(model.win3d, StsMessage.WARNING,
                            "Failed to find float data volume to read: switching to byte data .");
                    isInputFloats = false;
                    inputFormatBean.setValue(inputFormatByte);
                }
                inputFormatBean.setValue(inputFormatFloat);
            }
        }
        else if(!volume.setupRowByteBlocks())
        {
            new StsMessage(model.win3d, StsMessage.WARNING, "Failed to find byte data volume to read: cannot process.");
            return false;
        }

        isOutputFloats = (outputFormat == IEEE_FLOAT || outputFormat == IBM_FLOAT);
        outputTraceFloats = new float[nCroppedSlices];
        outputTraceBytes = new byte[nCroppedSlices];
        sampleSize = isOutputFloats ? 4 : 1;

        if(outputFormat == SCALED8 || outputFormat == BYTE)
            traceHeaderScale = unitsScale * volume.getByteScale();
        else
            traceHeaderScale = 1.0f;

        return true;
    }

    private boolean setGridParameters()
    {
        croppedBoundingBox = new StsCroppedBoundingBox(volume, exportVolume, false);
        nVolumeRows = volume.nRows;
        nVolumeCols = volume.nCols;
        nVolumeSlices = volume.nSlices;
        nCroppedSlices = croppedBoundingBox.nSlices;
        return true;
    }

    private boolean exportRow(StsSEGYFormat segyFmt, int row, float rowNum) throws Exception
    {
        byte[] inputBytes = null;
        float[] inputFloats = null;

        // Read in New Line and convert as necessary
        if(isInputFloats)
        {
            inputFloats = volume.readRowPlaneFloatData(row);
            if(inputFloats == null)
            {
                new StsMessage(model.win3d, StsMessage.WARNING, "Failed to read row plane data for volume(" + volume.getName() + ")");
                return false;
            }
        }
        else
        {
            inputBytes = volume.readRowPlaneByteData(row);
            if(inputBytes == null)
            {
                new StsMessage(model.win3d, StsMessage.WARNING, "Failed to read row plane data for volume(" + volume.getName() + ")");
                return false;
            }
        }
        int offset = 0;
        int colMin = croppedBoundingBox.colMin;
        int colMax = croppedBoundingBox.colMax;
        int colInc = croppedBoundingBox.colInc;
        float colNum = croppedBoundingBox.colNumMin;
        float colNumInc = croppedBoundingBox.colNumInc;
        float zInc = croppedBoundingBox.zInc;
        for(int col = colMin; col <= colMax; col += colInc, colNum += colNumInc, offset += nVolumeSlices)
        {
            double[] xy = volume.getAbsoluteXY(volume.getXCoor(col), volume.getYCoor(row));
            constructTraceHdr(segyFmt, rowNum, colNum, xy[0] * 100.0, xy[1] * 100.0, nCroppedSlices, zInc, 1.0f, row * nVolumeRows + col, -100.0f, traceHeaderScale);
            dos.write(traceHeader);

            if(isInputFloats)
            {
                convertInputFloatsToOutputFloats(inputFloats, outputTraceFloats, offset, nCroppedSlices, row, col);
                if(isOutputFloats)
                {
                    if(!writeOutputFloats(outputTraceFloats, nCroppedSlices))
                        return false;
                }
                else // outputBytes
                {
                    convertOutputFloatsToOutputBytes(outputTraceFloats, outputTraceBytes, nCroppedSlices);
                    if(!writeOutputBytes(outputTraceBytes))
                        return false;
                }
            }
            else // inputBytes
            {
                convertInputBytesToOutputBytes(inputBytes, outputTraceBytes, offset, nCroppedSlices);
                if(isOutputFloats)
                {
                    convertOutputBytesToOutputFloats(outputTraceBytes, outputTraceFloats);
                    if(!writeOutputFloats(outputTraceFloats, nCroppedSlices))
                        return false;
                }
                else // outputBytes
                {
                    if(!writeOutputBytes(outputTraceBytes))
                        return false;
                }
            }
        }
        return true;
    }

    static public void main(String[] args)
    {
        try
        {
            StsSeismicVolume volume = StsSeismicVolume.constructTestVolume(args[0]);
            StsSeismicExportPanel3d exportPanel = new StsSeismicExportPanel3d(StsModel.currentModel, volume, "3D Export", false);
            StsToolkit.createDialog("Test", exportPanel, true, 250, 500);
        }
        catch(Exception e)
        {
            StsException.outputWarningException(StsSeismicVolume3dGroupBox.class, "main", e);
        }
    }
}
