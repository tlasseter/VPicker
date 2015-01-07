package com.Sts.UI;

import com.Sts.DBTypes.*;
import com.Sts.IO.StsAsciiFile;
import com.Sts.IO.StsFile;
import com.Sts.MVC.StsModel;
import com.Sts.MVC.StsProject;
import com.Sts.Types.StsPoint;
import com.Sts.Types.StsRotatedGridBoundingBox;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.StsException;
import com.Sts.Utilities.StsParameters;
import com.Sts.Utilities.StsStringUtils;

import javax.swing.*;
import java.awt.*;
import java.io.File;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

public class StsFractureSetExportDialog extends JDialog
{
    StsModel model;
    StsFractureSet fractureSet = null;
	String exportName;
    boolean exportEdges = false;

	StsJPanel panel = StsJPanel.addInsets();
	StsGroupBox selectBox;
	StsStringFieldBean nameBean;
    StsBooleanFieldBean exportEdgesBean;

	StsJPanel buttonPanel = StsJPanel.addInsets();
	StsButton processButton = new StsButton("Process", "Export selected curves of selected fractureSet.", this, "process");
	StsButton cancelButton = new StsButton("Cancel", "Cancel this operation.", this, "cancel");

    public final static byte PROCESS = 0;
    public final static byte CANCELED = 1;
    byte mode = PROCESS;

    protected StsFractureSetExportDialog(StsModel model, Frame frame, String title, boolean modal, StsFractureSet fractureSet, String timeOrDepth)
    {
        super(frame, title, modal);
        this.model = model;
        this.fractureSet = fractureSet;
        exportName = fractureSet.getName();
        constructBeans();
        constructPanel();
    }

    protected void constructBeans()
    {
        selectBox = new StsGroupBox("Export fracture set");
        nameBean = new StsStringFieldBean(this, "exportName", "Exported fractureSet name");
        exportEdgesBean = new StsBooleanFieldBean(this, "exportEdges", "Export as Edges");
        exportEdgesBean.setValue(true);
        exportEdgesBean.setEditable(false);
    }

    static public boolean exportFractureSet(StsModel model, Frame frame, String title, boolean modal, StsFractureSet well, String timeOrDepth)
    {
        try
        {
            StsFractureSetExportDialog dialog = new StsFractureSetExportDialog(model, frame, title, modal, well, timeOrDepth);
            dialog.setSize(200, 400);
            dialog.pack();
            dialog.setVisible(true);
            dialog.exportFractureSet(timeOrDepth);
            return true;
        }
        catch(Exception e)
        {
            StsException.systemError(StsFractureSetExportDialog.class, "constructor");
            return false;
        }
    }

    protected void constructPanel()
    {

	    this.getContentPane().add(panel);
		this.setTitle("Well Export Parameters");
		panel.add(selectBox);
		panel.add(buttonPanel);
		selectBox.add(nameBean);
        selectBox.add(exportEdgesBean);

	    buttonPanel.addToRow(processButton);
		buttonPanel.addEndRow(cancelButton);
    }

	public void setExportName(String name)
    {
        exportName = StsStringUtils.cleanString(name);
    }
	public String getExportName() { return exportName; }

    public void setExportEdges(boolean exportEdges)
    {
        this.exportEdges = exportEdges;
    }
    public boolean getExportEdges() { return exportEdges; }
    public void process()
	{
		mode = PROCESS;
        setVisible(false);
	}

	public void cancel()
	{
		mode = CANCELED;
		setVisible(false);
 	}

    public byte getMode() { return mode; }

    public boolean exportFractureSet(String timeOrDepth)
    {
        StsAsciiFile asciiFile = null;

        if (getMode() == CANCELED)
        {
            return false;
        }

        try
        {
            StsProject project = model.getProject();
            String fractureSetPathname = project.getRootDirString() + "fractureSet.txt." + exportName;
            StsObjectRefList fractureList = fractureSet.getFractureList();
            int nFractures = fractureList.getSize();
            for(int n = 0; n < nFractures; n++)
            {
                StsFracture fracture = (StsFracture)fractureList.getElement(n);
                String fracturePathname = fractureSetPathname + "." + n;
                File file = new File(fracturePathname);
                if (file.exists())
                {
                    boolean overWrite = StsYesNoDialog.questionValue(model.win3d,
                        "File " + fracturePathname + " already exists. Do you wish to overwrite it?");
                    if (!overWrite) continue;
                }
                StsFile stsFile = StsFile.constructor(fracturePathname);
                asciiFile = new StsAsciiFile(stsFile);
                if(!asciiFile.openWrite())
                {
                    new StsMessage(model.win3d, StsMessage.WARNING, "Failed to open file for writing: " + fracturePathname);
                    return false;
                }

                if(exportEdges)
                {
                    StsObjectRefList sectionEdges = fracture.getSectionEdges();
                    int nSectionEdges = sectionEdges.getSize();
                    for(int e = 0; e < nSectionEdges; e++)
                    {
                        StsSectionEdge edge = (StsSectionEdge)sectionEdges.getElement(e);
                        StsPoint[] points = edge.getPoints();
                        int nPoints = points.length;
                        for(int p = 0; p < nPoints; p++)
                        {
                            StsPoint point = points[p];
                            int position = 2;
                            if(p == 0)
                                position = 1;
                            else if(p == nPoints-1)
                                position = 3;
                            double[] xy = project.getAbsoluteXYCoordinates(point);
                            double x = xy[0];
                            double y = xy[1];
                            float z = point.getZ();
                            float rowNum = getRowNumFromY(project, point.getY());
                            float colNum = getColNumFromX(project, point.getX());
                            asciiFile.writeLine(x + "\t" + y + "\t" + z + "\t" + rowNum + "\t" + colNum + "\t" + position);
                        }
                    }
                }
                else
                {

                }
                asciiFile.close();
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsFractureSet.export() failed for fractureSet " + getName(), e, StsException.WARNING);
            return false;
        }
        finally
        {
            if(asciiFile != null) asciiFile.close();
        }
    }

    private float getRowNumFromY(StsProject project, float y)
    {
        StsRotatedGridBoundingBox boundingBox = project.rotatedBoundingBox;
        if(boundingBox.rowNumMin == StsParameters.nullValue) return 0.0f;
        return boundingBox.getRowNumFromY(y);
    }

    private float getColNumFromX(StsProject project, float x)
    {
        StsRotatedGridBoundingBox boundingBox = project.rotatedBoundingBox;
        if(boundingBox.colNumMin == StsParameters.nullValue) return 0.0f;
        return boundingBox.getColNumFromX(x);
    }
}
