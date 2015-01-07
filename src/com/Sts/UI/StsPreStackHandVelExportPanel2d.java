package com.Sts.UI;

import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.Types.PreStack.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;
import com.Sts.Utilities.*;

import java.awt.*;
import java.text.*;

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

public class StsPreStackHandVelExportPanel2d extends StsPreStackHandVelExportPanel
{
	public StsPreStackHandVelExportPanel2d(StsModel model, StsPreStackVelocityModel volume)
	{
		super(model, volume);
	}

	public boolean exportHandVels(StsProgressPanel progressPanel)
	{
		DecimalFormat labelFormat = new DecimalFormat("###0.0");

		StsAsciiFile asciiFile = null;
		String tempStg = new String();

		StsVelocityProfile[] velocityProfiles = (StsVelocityProfile[])volume.getVelocityProfiles().getCastList();
		if (velocityProfiles == null || velocityProfiles.length <= 0)
		{
			new StsMessage(model.win3d, StsMessage.WARNING, "No velocity profiles exist for this volume.");
			return false;
		}
		
		pathnames = new String[lineSet.lines.length];

		// Export the velocity profiles
		try
		{
			int nVelocityProfiles;
			int lastRow = -1;

			// Output of actual hand picked velocity profiles
//			if (getAllOrSomeProfiles() != INTERPOLATED_PROFILES)
			{
				nVelocityProfiles = velocityProfiles.length;
				progressPanel.initialize(nVelocityProfiles);
				for (int i = 0; i < nVelocityProfiles; i++)
				{
					StsVelocityProfile profile = velocityProfiles[i];
					//this problem should be fixed now (june 29, 2009), but still sort to make sure handvels will be OK in FOCUS
					profile.sort();
					profile.unique();
					profile.unique(); //unique twice in case there are 3 picks at same time.
					if (asciiFile == null || lastRow != profile.row)
					{
						if (asciiFile != null) {
							asciiFile.writeLine(getFileFooter());
							asciiFile.close();
						}
						lastRow = profile.row;
                        pathnames[lastRow] = getOutputPathname(volume.lineSet.lines[lastRow]);
						asciiFile = setupExport(pathnames[lastRow], volume.lineSet.lines[lastRow]);
						if (asciiFile == null)
							return false;
					}
                    int cdpNum = lineSet.getCDP(profile.row, profile.col);
//                    float profileRowNum = lineSet.getRowNumFromRow(profile.row);
//					float profileColNum = lineSet.getColNumFromCol(profile.col);
//					if ((profileRowNum >= rowNumMin && profileRowNum <= rowNumMax) && (profileColNum >= colNumMin && profileColNum <= colNumMax))
					{
						// Export CDP number
//						int cdpNum = (int)((volume.getRowNumFromRow(profile.row) * 10000.0f) + volume.getColNumFromCol(profile.col));
						tempStg = "HANDVEL " + cdpNum;
						tempStg = StsStringUtils.padClipString(tempStg, 80);
						asciiFile.writeLine(tempStg);
						// Export Time/Velocity Pairs
						tempStg = "";
						StsPoint[] semblancePoints = profile.getProfilePoints();
						int nSemblancePoints = semblancePoints.length;
						for (int j = 0; j < nSemblancePoints; j++)
						{
							float[] flts = semblancePoints[j].v;
							String location = StsStringUtils.padClipString(Integer.toString((int)flts[1]), 8);
							String velocity = StsStringUtils.padClipString(labelFormat.format(flts[0] * unitsScale), 8);

							tempStg = tempStg + location + velocity;
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
        /*
            else
			{
				FloatBuffer floatBuffer = null;
				float[] traceFloats = null;
				int nCroppedSlices = volume.nCroppedSlices;

				nVelocityProfiles = (int)Math.abs((rowNumMax - rowNumMin) / rowNumInc) * (int)Math.abs((colNumMax - colNumMin) / colNumInc);
				progressPanel.initialize(nVelocityProfiles);

				traceFloats = new float[nCroppedSlices];
				for (float rowNum = rowNumMin; rowNum < rowNumMax; rowNum += rowNumInc)
				{
					int row = volume.getRowFromRowNum(rowNum);

					if (asciiFile != null)
						asciiFile.close();
					String pathname = project.getRootDirString() + volume.getName() + "-line" + row + ".handvel";
					asciiFile = setupExport(pathname);
					if (asciiFile == null)
						return false;

					floatBuffer = ((StsPreStackVelocityModel2d)volume).computeFloatBufferVelocityPlane(row);
					int maxCol = ((StsPreStackLineSet2d)this.lineSet).lines[row].nCols;
					for (float colNum = colNumMin; colNum < maxCol; colNum += colNumInc)
					{
						int cdpNum = (int)((rowNum * 10000.0f) + colNum);
						tempStg = "HANDVEL " + cdpNum;
						tempStg = StsStringUtils.padClipString(tempStg, 80);
						asciiFile.writeLine(tempStg);

						int col = volume.getColFromColNum(colNum);
						floatBuffer.position(col * nCroppedSlices);
						floatBuffer.get(traceFloats);
						int cnt = 0;
						tempStg = "";
						for (float n = zMin; n < nCroppedSlices; n += zInc)
						{
							int traceNum = (int)(n / volume.getZInc());
							float traceFloat = traceFloats[traceNum];
							String location = StsStringUtils.padClipString(Integer.toString((int)n), 8);
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
						progressPanel.incrementProgress();
					}
				}
			}
        */
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
					
				} catch (Exception e1)
				{
					e1.printStackTrace();
				}
			}
		}

		return true;
	}
	
	@Override
	protected String getOutputPathname(StsPreStackLine line)
    {
	    StsProject project = model.getProject();
	    String lineName = line.getName();
	    if( lineName.length() > 0)
            return project.getRootDirString() + getExportName() + "_" + lineName + "_handvel.dat"; //Focus wants filenames with only one "." and end in ".dat" SWC 6/3/09
        else
            return project.getRootDirString() + getExportName() + "_line" + line.getRowNumMin() + "_handvel.dat";
    }

	@Override
    protected String getFileFooter()
    {
        return "__End Handvels__";
    }

    protected void assembleGroupBoxes() throws Exception
	{
		exportParametersGroupBox.addToRow(exportNameBean);
		exportParametersGroupBox.addEndRow(nullValueBean);
		exportParametersGroupBox.addToRow(velocityTypeBean);
		exportParametersGroupBox.addEndRow(exportDomainBean);
		exportParametersGroupBox.addToRow(templateBean, 2, 1.0);

		constructVelocityUnitsGroupBox();
	}

	protected void assemblePanel()
	{
		titleLabel.setFont(new java.awt.Font("Serif", 1, 14));
		gbc.fill = GridBagConstraints.HORIZONTAL;
		add(titleLabel);
//		add(profileButtons.getButtonPanel());
		add(exportParametersGroupBox);
		add(unitsGroupBox);
		gbc.fill = GridBagConstraints.NONE;
		add(textHdrEditBtn);
		gbc.fill = GridBagConstraints.HORIZONTAL;
		//add(progressPanel);
	}
}
