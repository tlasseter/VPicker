package com.Sts.DBTypes;

import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.io.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2007</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */

/** for 2d we have a set of lines in a volume (lineVolume); each line is exported separately and is treated as a volume
 *  (StsSeismicBoundingBox) in order to reuse logic in StsSeismicExportPanel.  So for each line in the lineVolume, we
 *  set volume equal to that line and all else follows.
 */
public class StsSeismicExportPanel2d extends StsSeismicExportPanel
{
    StsSeismicBoundingBox lineVolume;
    /** col range for exported volume */
	transient int colStart, colEnd, colInc;

    public StsSeismicExportPanel2d()
	{
		super();
	}

	public StsSeismicExportPanel2d(StsModel model, StsSeismicBoundingBox volume, String title, boolean isAVelocityExport)
	{
		super(model, volume, title, isAVelocityExport);
        lineVolume = volume;
    }

	protected void buildBeans()
	{
		nRows = volume.getNRows();
		super.buildBeans();
	}

    protected StsGroupBox getInputGeometryGroupBox()
    {
        return new StsSeismicVolume2dGroupBox(volume, "Selected Volume");
    }

    protected StsGroupBox getOutputGeometryGroupBox()
    {
        return new StsSeismicVolume2dGroupBox(exportVolume, "Export Volume", volume);
    }

	protected void doExportSeismic(StsProgressPanel progressPanel)
	{
		if (!setupExport())
        {
            new StsMessage(model.win3d,  StsMessage.WARNING, "Failed to setup for export.");
            return;
        }

        StsSEGYFormat segyFmt = StsSEGYFormat.constructor(model, StsSEGYFormat.POSTSTACK);

        byte currentType = velocityType;
        // type is not saved in segy, so we need to find a place for it in the segy file headers
        lineVolume.setType(StsParameters.SAMPLE_TYPE_VEL_RMS);

        progressPanel.setMaximum(nRows);
        for (int row = 0; row < nRows && ! canceled; row++)
        {
            String filename = exportName + ".line." + row + ".sgy";
            if (lineVolume.hasExportableData(row))
            {
                try
                {
                    exportRow(segyFmt, row, filename);
                }
                catch (Exception e)
                {
                    progressPanel.appendLine("Export failed for " + filename);
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
                    catch (Exception e)
                    {
                        progressPanel.appendLine("Export failed for " + filename);
                        StsException.outputException("StsSeismicExportPanel2d.doExportSeismic() failed.", e, StsException.WARNING);
                        return;
                    }
                }
                progressPanel.appendLine("Exported " + filename);
            }
            else
            {
                progressPanel.appendLine("No data to export for " + filename);
            }
            progressPanel.setValue(row+1);
        }

        // type is not saved in segy, so we need to find a place for it in the segy file headers
        lineVolume.setType(currentType);
	}

	protected boolean checkExportSeismic()
	{
        boolean checked = false;
        boolean overwrite = false;
        for(int row = 0; row < nRows; row++)
		{
			String filename = model.getProject().getRootDirString() + exportName + ".line." + row + ".sgy";
			File file = new File(filename);
			if(file.exists())
            {
               if(!checked)
			   {
                    overwrite = StsYesNoDialog.questionValue(model.win3d, "One or more 2d lines exist. Do you wish to overwrite them?");
				    if (!overwrite) return false;
                    checked = true;
               }
               if(overwrite) file.delete();
		    }
        }
        return true;
	}

	protected boolean setupExport()
	{
        isInputFloats = lineVolume.isDataFloat;
		if (isInputFloats)
		{
			if (!lineVolume.setupReadRowFloatBlocks())
			{
				if (!lineVolume.setupRowByteBlocks())
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
			}
            inputFormatBean.setValue(inputFormatFloat);
		}
		else if (!lineVolume.setupRowByteBlocks())
		{
			new StsMessage(model.win3d, StsMessage.WARNING, "Failed to find byte data volume to read: cannot process.");
			return false;
		}

		isOutputFloats = (outputFormat == IEEE_FLOAT || outputFormat == IBM_FLOAT);
//        outputTraceFloats = new float[nCroppedSlices];
//        outputTraceBytes = new byte[nCroppedSlices];
		sampleSize = isOutputFloats ? 4 : 1;

		if (outputFormat == SCALED8 || outputFormat == BYTE)
			traceHeaderScale = unitsScale * lineVolume.getByteScale();
		else
			traceHeaderScale = 1.0f;

		return true;
	}

    private void exportRow(StsSEGYFormat segyFmt, int row, String filename) throws Exception
	{
        byte[] inputBytes = null;
		float[] inputFloats = null;

        String pathname = model.getProject().getRootDirString() + filename;
        setupFile(segyFmt, pathname);

        setGridParameters(row);

        constructTextHeader();
		constructBinaryHdr();

        writeTextHdr(segyFmt);
        writeBinaryHeader();

//		byte[] outputBytes = null;
//		float[] outputFloats = null;

        // Read in New Line and convert as neccessary
		if (isInputFloats)
			inputFloats = lineVolume.readRowPlaneFloatData(row);
		else
			inputBytes = lineVolume.readRowPlaneByteData(row);

        if(isOutputFloats)
            outputTraceFloats = new float[nCroppedSlices];
        else
            outputTraceBytes = new byte[nCroppedSlices];
        
        double[][] cdps = lineVolume.getLineXYCDPs(row);
        int offset = 0;
        int nSlices = croppedBoundingBox.nSlices;
		for(int col = 0; col < nCols; col++, offset += nVolumeSlices)
		{
			constructTraceHdr(segyFmt, row, col, cdps[col][0] * 100.0, cdps[col][1] * 100.0, nSlices, exportVolume.zInc, 1.0f,
									cdps[col][2], -100.0f, traceHeaderScale);
			dos.write(traceHeader);
			if (isInputFloats)
			{
                convertInputFloatsToOutputFloats(inputFloats, outputTraceFloats, offset, nSlices, row, col);
                if (isOutputFloats)
				{
                    if (!writeOutputFloats(outputTraceFloats, nSlices))
						return;
				}
				else
				{

                    convertOutputFloatsToOutputBytes(outputTraceFloats, outputTraceBytes, nSlices);
                    if (!writeOutputBytes(outputTraceBytes))
						return;
				}
			}
			else // inputBytes
			{
                convertInputBytesToOutputBytes(inputBytes, outputTraceBytes, offset, nSlices);
                if (isOutputFloats)
				{
                    convertOutputBytesToOutputFloats(outputTraceBytes, outputTraceFloats);
                    if (!writeOutputFloats(outputTraceFloats, nSlices))
						return;
				}
				else
				{
					if (!writeOutputBytes(outputTraceBytes))
						return;
				}
			}
		}
	}

    //TODO adjust code to reflect user changes from GUI; for now we are ignoring them and exporting all samples for the entire volume
    private boolean setGridParameters(int row)
    {
        volume = lineVolume.getLineBoundingBox(row);
		nVolumeCols = lineVolume.getColumnsInRow(row);
        nCols = nVolumeCols;
        exportVolume.zInc = lineVolume.getZInc(row);
        sliceStart = 0;
        sliceInc = 1;
        nVolumeSlices = lineVolume.getNSlices(row);
        nCroppedSlices = nVolumeSlices;
        return true;
    }
}
