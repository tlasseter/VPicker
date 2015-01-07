package com.Sts.UI;

import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.Types.PreStack.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Progress.*;
import com.Sts.Utilities.*;

import java.io.*;
import java.nio.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2007</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version 1.0
 */

public class StsPreStackHandVelExportPanel3d extends StsPreStackHandVelExportPanel
{
    public static final String CDP = "CDP";
    public static final String CDPLBL = "CDPLBL";
    public static final String[] SortOrders = new String[]{ CDPLBL, CDP };
    public String cdpSortOrder = CDPLBL;
    
    protected StsComboBoxFieldBean cdpSortOrderBean;
    
	public StsPreStackHandVelExportPanel3d(StsModel model, StsPreStackVelocityModel volume)
	{
		super(model, volume);
	}

	public boolean exportHandVels(StsProgressPanel progressPanel)
	{
		StsDecimalFormat labelFormat = new StsDecimalFormat(5);

		StsAsciiFile asciiFile = null;
		String tempStg;

		StsVelocityProfile[] velocityProfiles = (StsVelocityProfile[])volume.getVelocityProfiles().getCastList();
		if (velocityProfiles == null || velocityProfiles.length <= 0)
		{
			new StsMessage(model.win3d, StsMessage.WARNING, "No velocity profiles exist for this volume.");
			return false;
		}

		// Export the velocity profiles
		pathnames = new String[] { getOutputPathname(lineSet.lines[0]) };
		
		try
		{
			asciiFile = setupExport(pathnames[0], lineSet.lines[0]);
			if (asciiFile == null)
				return false;

			int nVelocityProfiles;

			//this test necessary for "backwards" 3D with negative inline increment
			float minTmp = rowNumMin;
			rowNumMin = Math.min(rowNumMin, rowNumMax);
			rowNumMax = Math.max(minTmp, rowNumMax);
			minTmp = colNumMin;
			colNumMin = Math.min(colNumMin, colNumMax);
			colNumMax = Math.max(minTmp, colNumMax);
			
			// Output of actual hand picked velocity profiles
			if (getAllOrSomeProfiles() != INTERPOLATED_PROFILES)
			{
				nVelocityProfiles = velocityProfiles.length;
				progressPanel.initialize(nVelocityProfiles);
				for (int i = 0; i < nVelocityProfiles; i++)
				{
					if (canceled)
						return false;
					StsVelocityProfile profile = velocityProfiles[i];
					//this problem should be fixed now (june 29, 2009), but still sort to make sure handvels will be OK in FOCUS
					profile.sort();
					profile.unique();
					profile.unique(); //unique twice in case there are 3 picks at same time.
					float profileRowNum = volume.getRowNumFromRow(profile.row);
					float profileColNum = volume.getColNumFromCol(profile.col);
					if ((profileRowNum >= rowNumMin && profileRowNum <= rowNumMax) && (profileColNum >= colNumMin && profileColNum <= colNumMax))
					{
						// Export CDP number
						tempStg = "HANDVEL " + getCDP(profile.row, profile.col);
						tempStg = StsStringUtils.padClipString(tempStg, 80);
						asciiFile.writeLine(tempStg);
						// Export Time/Velocity Pairs
						tempStg = "";
						StsPoint[] semblancePoints = profile.getProfilePoints();
						int nSemblancePoints = semblancePoints.length;
						for (int j = 0; j < nSemblancePoints; j++)
						{
							float[] flts = semblancePoints[j].v;
							String velocity = StsStringUtils.padClipString(labelFormat.formatValue(flts[0] * unitsScale), 8);
                            String timeDepth = StsStringUtils.padClipString(Integer.toString((int)flts[1]), 8);

                            tempStg = tempStg + timeDepth + velocity;
							if (((j + 1) % 4 == 0) && (j != 0))
							{
								tempStg = StsStringUtils.padClipString(tempStg, 80);
								asciiFile.writeLine(tempStg);
								tempStg = "";
							}
						}
						if (tempStg != "")
						{
							tempStg = StsStringUtils.padClipString(tempStg, 80);
							asciiFile.writeLine(tempStg);
						}
					}
					progressPanel.setValue(i+1);
				}
			}
			// Output of interpolated profiles at user specified intervals
			else
			{
				FloatBuffer floatBuffer = null;
				float[] traceFloats = null;
				int nSlices = volume.nSlices;

				nVelocityProfiles = (int)Math.abs((rowNumMax - rowNumMin) / rowNumInc) * (int)Math.abs((colNumMax - colNumMin) / colNumInc);
				progressPanel.initialize(nVelocityProfiles);
                int n = 0;
				traceFloats = new float[nSlices];
				for (float rowNum = rowNumMin; rowNum < rowNumMax; rowNum += rowNumInc)
				{
					int row = volume.getRowFromRowNum(rowNum);
					floatBuffer = ((StsPreStackVelocityModel3d)volume).computeFloatBufferVelocityPlane(StsParameters.YDIR, row);
					for (float colNum = colNumMin; colNum < colNumMax; colNum += colNumInc)
					{
						if (canceled)
							return false;
						int col = volume.getColFromColNum(colNum);
						tempStg = "HANDVEL " + getCDP(row, col);
						tempStg = StsStringUtils.padClipString(tempStg, 80);
						asciiFile.writeLine(tempStg);

						floatBuffer.position(col * nSlices);
						floatBuffer.get(traceFloats);
						int cnt = 0;
						tempStg = "";
						for (float z = zMin; z < zMax; z += zInc)
						{
							int traceNum = (int)(n / volume.getZInc());
							float traceFloat = traceFloats[traceNum];
                            int zi = Math.round(z);
                            String location = StsStringUtils.padClipString(Integer.toString(zi), 8);
							String velocity = StsStringUtils.padClipString(labelFormat.format(traceFloat * unitsScale), 8);
							tempStg = tempStg + location + velocity;
							if (((cnt + 1) % 4 == 0) && (cnt != 0))
							{
								tempStg = StsStringUtils.padClipString(tempStg, 80);
								asciiFile.writeLine(tempStg);
								tempStg = "";
							}
							cnt++;
						}
						if (tempStg != "")
						{
							tempStg = StsStringUtils.padClipString(tempStg, 80);
							asciiFile.writeLine(tempStg);
							tempStg = "";
						}
						progressPanel.setValue(++n);
					}
				}
			}
		}
		catch (Exception e)
		{
			StsException.outputException("StsPreStackVelocityModel.export() failed for volume " + volume.getName(), e, StsException.WARNING);
			return false;
		}
		finally
		{
		    if (asciiFile != null)
		    {
		        try
                {
                    asciiFile.writeLine(getFileFooter());
                    asciiFile.close();
                } catch (IOException e)
                {
                    e.printStackTrace();
                }
		    }
		}
		return true;
	}

    protected String getCDP(int row, int col)
    {
        int rowNum = (int) volume.getRowNumFromRow(row);
        int colNum = (int) volume.getColNumFromCol(col);
        if (cdpSortOrder.equals(CDP)) return "" + lineSet.getCDPFromRowNumColNum(rowNum, colNum);
        int cdpNum = rowNum * 10000 + colNum;
        return "" + cdpNum;
    }

    protected String getFileFooter()
    {
        return "__End Handvels__";
    }

    @Override
    protected void assembleGroupBoxes() throws Exception
    {
        exportParametersGroupBox.addToRow(exportNameBean);
        exportParametersGroupBox.addEndRow(nullValueBean);
        exportParametersGroupBox.addToRow(velocityTypeBean);
        exportParametersGroupBox.addEndRow(exportDomainBean);
        exportParametersGroupBox.addToRow(templateBean);
        exportParametersGroupBox.addEndRow(cdpSortOrderBean);

        subVolumeGroupBox.gbc.gridx = 1; // leave empty cell
        subVolumeGroupBox.addToRow(minLabel);
        subVolumeGroupBox.addToRow(maxLabel);
        subVolumeGroupBox.addEndRow(intervalLabel);

        subVolumeGroupBox.addToRow(rowNumMinBean);
        subVolumeGroupBox.addToRow(rowNumMaxBean);
        subVolumeGroupBox.addEndRow(rowNumIncBean);

        subVolumeGroupBox.addToRow(colNumMinBean);
        subVolumeGroupBox.addToRow(colNumMaxBean);
        subVolumeGroupBox.addEndRow(colNumIncBean);

        subVolumeGroupBox.addToRow(zMinBean);
        subVolumeGroupBox.addToRow(zMaxBean);
        subVolumeGroupBox.addEndRow(zIncBean);

        constructVelocityUnitsGroupBox();
        
    }
    
    @Override
    protected void buildBeans()
    {
        super.buildBeans();
        cdpSortOrderBean = new StsComboBoxFieldBean(this, "cdpSortOrder", "Sort Order", SortOrders );
    }
    
    public String getCdpSortOrder()
    {
        return cdpSortOrder.toString();
    }

    public void setCdpSortOrder(String cdpSortOrder)
    {
        this.cdpSortOrder = cdpSortOrder;
    }

    @Override
    protected String getOutputPathname(StsPreStackLine line)
    {
        StsProject project = model.getProject();
        return project.getRootDirString() + getExportName() + "_" + project.getName() + "_handvel.dat";
    }
	
}
