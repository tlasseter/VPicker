package com.Sts.UI;

import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.xml.bind.JAXBContext;
import javax.xml.bind.Marshaller;
import java.awt.*;
import java.io.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

public class StsWellExportDialog extends JDialog
{
    StsModel model;
    StsWell well = null;
    boolean logsData = false;
    boolean deviationData = true;
    boolean seisAttsData = false;
    float minMDepth;
	float maxMDepth;
    float sampleRate;
    boolean resample = false;
	String exportName;

	StsJPanel panel = StsJPanel.addInsets();
	StsGroupBox selectBox;
	StsStringFieldBean nameBean;
	StsBooleanFieldBean logsDataBean;
	StsBooleanFieldBean deviationDataBean;
	StsBooleanFieldBean seisAttsDataBean;

    StsGroupBox rangeBox;
	StsFloatFieldBean minMDepthBean;
	StsFloatFieldBean maxMDepthBean;
    StsBooleanFieldBean resampleBean;
    StsFloatFieldBean sampleRateBean;

	StsJPanel buttonPanel = StsJPanel.addInsets();
	StsButton processButton = new StsButton("Process", "Export selected curves of selected well.", this, "process");
	StsButton cancelButton = new StsButton("Cancel", "Cancel this operation.", this, "cancel");

    public final static byte PROCESS = 0;
    public final static byte CANCELED = 1;
    byte mode = PROCESS;

    protected StsWellExportDialog(StsModel model, Frame frame, String title, boolean modal, StsWell well, String timeOrDepth)
    {
        super(frame, title, modal);
        this.model = model;
        this.well = well;
        exportName = well.getName();
        constructBeans();
        constructPanel();
    }

    protected void constructBeans()
    {
        selectBox = new StsGroupBox("Export curves");
        nameBean = new StsStringFieldBean(this, "exportName", "Exported well name");
        logsDataBean = new StsBooleanFieldBean(this, "logsData", "Log Data");
        deviationDataBean = new StsBooleanFieldBean(this, "deviationData", "Dev Curve");
        seisAttsDataBean = new StsBooleanFieldBean(this, "seisAttsData", "Seismic Attributes");

        rangeBox = new StsGroupBox("Export mdepth range");
        minMDepthBean = new StsFloatFieldBean(this, "minMDepth", "Min MDepth");
        maxMDepthBean = new StsFloatFieldBean(this, "maxMDepth", "Max MDepth");

        resampleBean = new StsBooleanFieldBean(this, "resample", "Resample Output");
        sampleRateBean = new StsFloatFieldBean(this, "sampleRate", "Samples Per Foot");
    }

    static public boolean exportWell(StsModel model, Frame frame, String title, boolean modal, StsWell well, String timeOrDepth)
    {
        try
        {
            StsWellExportDialog dialog = new StsWellExportDialog(model, frame, title, modal, well, timeOrDepth);
            dialog.setSize(200, 400);
            dialog.pack();
            dialog.setVisible(true);
            dialog.exportWell(timeOrDepth);
            return true;
        }
        catch(Exception e)
        {
            StsException.systemError(StsWellExportDialog.class, "constructor");
            return false;
        }
    }

    protected void constructPanel()
    {

	    this.getContentPane().add(panel);
		this.setTitle("Well Export Parameters");
		panel.add(selectBox);
        panel.add(rangeBox);
		panel.add(buttonPanel);

		selectBox.add(nameBean);
		selectBox.addToRow(logsDataBean);
		selectBox.addToRow(deviationDataBean);
		selectBox.addEndRow(seisAttsDataBean);

        rangeBox.addToRow(minMDepthBean);
		rangeBox.addEndRow(maxMDepthBean);

        rangeBox.addToRow(resampleBean);
        rangeBox.addEndRow(sampleRateBean);

	    buttonPanel.addToRow(processButton);
		buttonPanel.addEndRow(cancelButton);

	    minMDepthBean.setValue(0.0f);
        minMDepth = 0.0f;
        maxMDepth = well.getLastPoint().getM();
		maxMDepthBean.setValue(maxMDepth);
        sampleRate = (maxMDepth = minMDepth)/well.getLineVertices().getSize();
        sampleRateBean.setValue(sampleRate);
        sampleRateBean.setRange(0.5, sampleRate);
        sampleRateBean.setEditable(false);
    }

	public void setExportName(String name)
    {
        exportName = StsStringUtils.cleanString(name);
    }
	public String getExportName() { return exportName; }

	public void setResample(boolean resample)
    {
        this.resample = resample;
        sampleRateBean.setEditable(resample);
    }
	public boolean getResample() { return resample; }

	public void setLogsData(boolean logsData) { this.logsData = logsData; }
	public boolean getLogsData() { return logsData; }

	public void setDeviationData(boolean deviationData) { this.deviationData = deviationData; }
	public boolean getDeviationData() { return deviationData; }

	public void setSeisAttsData(boolean seisAttsData) { this.seisAttsData = seisAttsData; }
	public boolean getSeisAttsData() { return seisAttsData; }

	public void setMinMDepth(float  minMDepth) { this.minMDepth = minMDepth; }
	public float getMinMDepth() { return this.minMDepth; }

	public void setMaxMDepth(float  maxMDepth) { this.maxMDepth = maxMDepth; }
	public float getMaxMDepth() { return this.maxMDepth; }

	public void setSampleRate(float sampleRate) { this.sampleRate = sampleRate; }
	public float getSampleRate() { return this.sampleRate; }

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

    public boolean exportWellXML(String timeOrDepth)
    {
        if (getMode() == CANCELED) return false;

        try
        {
            StsProject project = model.getProject();
            String pathname = project.getRootDirString() + "well-dev.xml." + exportName;
            File file = new File(pathname);
            if (file.exists())
            {
                boolean overWrite = StsYesNoDialog.questionValue(model.win3d,
                    "File " + pathname + " already exists. Do you wish to overwrite it?");
                if (!overWrite) return false;
            }
            JAXBContext context = JAXBContext.newInstance(StsWell.class, StsObjectRefList.class, StsLogCurve.class,
                    StsObjectList.class, StsSection.class, StsColor.class, StsSurfaceVertex.class, StsGridSectionPoint.class, StsBlock.class, StsSurface.class);
            Marshaller m = context.createMarshaller();
            m.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);
            m.marshal(this, file);
            // StsToolkit.writeObjectXML(well, pathname);
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsWell.export() failed for well " + getName(), e, StsException.WARNING);
            return false;
        }
    }

    public boolean exportWell(String timeOrDepth)
    {
        StsAsciiFile asciiFile = null;
        boolean exportLogData = getLogsData();
        boolean exportDeviationData = getDeviationData();
        boolean exportSeisAttData = getSeisAttsData();
		float minMDepth = getMinMDepth();
		float maxMDepth = getMaxMDepth();
        if (getMode() == CANCELED)
        {
            return false;
        }

        boolean writeTime = false;
        boolean writeDepth = false;
        if (timeOrDepth == StsParameters.TD_TIME_DEPTH_STRING)
        {
            writeDepth = true;
            writeTime = true;
        }
        else if (timeOrDepth == StsParameters.TD_DEPTH_STRING)
        {
            writeDepth = true;
        }
        else
        {
            writeTime = true;
        }

        try
        {
            StsProject project = model.getProject();
            String directory = project.getRootDirString();
            String pathname = project.getRootDirString() + "well-dev.txt." + exportName;
            File file = new File(pathname);
            if (file.exists())
            {
                boolean overWrite = StsYesNoDialog.questionValue(model.win3d,
                    "File " + pathname + " already exists. Do you wish to overwrite it?");
                if (!overWrite)
                {
                    return false;
                }
            }
            StsFile stsFile = StsFile.constructor(pathname);
            asciiFile = new StsAsciiFile(stsFile);
            if(!asciiFile.openWrite())
			{
				new StsMessage(model.win3d, StsMessage.WARNING, "Failed to open file for writing: " + pathname);
				return false;
			}
            asciiFile.writeLine("WELLNAME");
            asciiFile.writeLine(exportName);
            asciiFile.writeLine("ORIGIN XY");
            if(well instanceof StsWellPlan)
                asciiFile.writeLine(project.getXOrigin() + " " + project.getYOrigin());
            else
                asciiFile.writeLine(well.getXOrigin() + " " + well.getYOrigin());

            asciiFile.writeLine("CURVE");
            asciiFile.writeLine("X");
            asciiFile.writeLine("Y");
            if (writeDepth)
            {
                asciiFile.writeLine("DEPTH");
                asciiFile.writeLine("MDEPTH");
            }
            if (writeTime)
            {
                asciiFile.writeLine("TIME");
            }
//            int nVertices = lineVertices.getSize();
            // Header for Log Data
            StsObjectRefList curves = null;
            StsLogCurve log = null;
            if(exportLogData)
            {
                curves = well.getLogCurves();
                if(curves.getSize() > 0)
                {
                    for(int i=0; i<curves.getSize(); i++)
                    {
                        asciiFile.writeLine(curves.getElement(i).getName());
//                        if(((StsLogCurve)curves.getElement(i)).getValuesFloatVector().getSize() > nVertices)
//                        {
//                            log = (StsLogCurve)curves.getElement(i);
//                            nVertices = log.getValuesFloatVector().getSize();
//                        }
                    }
                }
            }
            // Header for Seismic Data
            StsSeismicVolume[] vols = null;
            if(exportSeisAttData)
            {
                vols = (StsSeismicVolume[])model.getCastObjectList(StsSeismicVolume.class);
                for(int i=0; i<vols.length; i++)
                    asciiFile.writeLine(vols[i].getName());
            }
            asciiFile.writeLine("VALUE");

//            String valLine = null;
			StsPoint[] points = well.getExportPoints();
			int nPoints = points.length;

            float exportMinDepth = minMDepth;
            if(well instanceof StsWellPlan)
            {
                StsWell drillingWell = ((StsWellPlan)well).getDrillingWell();
                if(drillingWell != null)
                    exportMinDepth = drillingWell.getMaxMDepth();
            }

            StsPoint minPoint = well.getPointAtMDepth(exportMinDepth, points, true);
            outputPoint(minPoint, writeTime, writeDepth, exportLogData, curves, asciiFile);
            if(resample)
                points = resamplePoints(points);
            for(int n = 0; n < points.length; n++)
            {
                float mdepth = points[n].getM();
                if(mdepth <= exportMinDepth) continue;
                if(mdepth >= maxMDepth) break;
                outputPoint(points[n], writeTime, writeDepth, exportLogData, curves, asciiFile);
            }
            StsPoint maxPoint = well.getPointAtMDepth(maxMDepth, points, true);
            outputPoint(maxPoint, writeTime, writeDepth, exportLogData, curves, asciiFile);

			if(exportSeisAttData)
			{
				;
			}
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsWell.export() failed for well " + getName(), e, StsException.WARNING);
            return false;
        }
        finally
        {
            if(asciiFile != null) asciiFile.close();
        }
    }

    private StsPoint[] resamplePoints(StsPoint[] points)
    {
        float minMDepth = getMinMDepth();
		float maxMDepth = getMaxMDepth();
        int numberSamples = (int)((maxMDepth - minMDepth) * getSampleRate());
        float inc = (int)((maxMDepth - minMDepth)/numberSamples);
        StsPoint[] newPts = new StsPoint[numberSamples];
        for(int i=0; i<numberSamples; i++)
        {
            newPts[i] = well.getPointAtMDepth(minMDepth + ((i+1)*inc), points, true);
        }
        return newPts;
    }
	private void outputPoint(StsPoint point, boolean writeTime, boolean writeDepth, boolean exportLogData, StsObjectRefList curves, StsAsciiFile asciiFile)
	{
		String valLine = null;

		float m = point.getM();
		float x = point.getX();
		float y = point.getY();
		float[] unrotatedXY = model.getProject().getUnrotatedRelativeXYFromRotatedXY(x, y);
		x = unrotatedXY[0];
		y = unrotatedXY[1];
		float z = point.getZ();
		float t = point.getT();
		if (writeTime && writeDepth)
			 valLine = new String(x + " " + y + " " + z + " " + m + " " + t);
		 else if (writeDepth)
			 valLine = new String(x + " " + y + z + " " + m);
		 else if (writeTime)
			 valLine = new String(x + " " + y + " " + t);

		if(exportLogData)
		{
			for(int i=0; i < well.getNLogCurves(); i++)
				valLine = valLine + " " + ((StsLogCurve)curves.getElement(i)).getInterpolatedValue(z);
		}
		try
		{
			asciiFile.writeLine(valLine);
		}
		catch(Exception e)
		{
		}
	}

    static public void main(String[] args)
	{
        StsModel model = StsModel.constructor();
        StsWell well = new StsWell();
		StsWellExportDialog exportDialog = new StsWellExportDialog(model, null, "Well Export Utility", true, well, "");
		 exportDialog.setVisible(true);
		 boolean exportLogData = exportDialog.getLogsData();
		 boolean exportDeviationData = exportDialog.getDeviationData();
		 boolean exportSeisAttData = exportDialog.getSeisAttsData();
		 float minMDepth = exportDialog.getMinMDepth();
		 float maxMDepth = exportDialog.getMaxMDepth();
	}
}
