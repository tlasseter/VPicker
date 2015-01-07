package com.Sts.Actions.Wizards.LogToVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSelectLogPanel extends StsJPanel
{
    private StsLogToVolumeWizard wizard;
    private StsSelectLog wizardStep;

    private StsModel model = null;

    StsLogCurve selectedLog = null;
    StsComboBoxFieldBean logBean = new StsComboBoxFieldBean();
    StsGroupBox logBox = new StsGroupBox("Select Log");
    StsLogCurve[] availableLogCurves = null;

    public static byte INSTANT = 0;
    public static byte AVERAGE = 1;
    public static byte[] computeMethods = new byte[] {INSTANT, AVERAGE};
    public static String[] computeMethodStrings = new String[] {"Instant", "Average"};
    public byte computeMethod = INSTANT;

    StsBooleanFieldBean resampleBean;
    StsIntFieldBean intervalBean;
    StsComboBoxFieldBean computeMethodBean;
    StsGroupBox computeBox = new StsGroupBox("ReSample Log");

    public StsSelectLogPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
    	this.wizard = (StsLogToVolumeWizard)wizard;
    	this.wizardStep = (StsSelectLog)wizardStep;
    	model  = wizard.getModel();
        if(!buildPanel())
        	wizard.cancel();
    }

    public boolean buildPanel()
    {
        gbc.fill = gbc.HORIZONTAL;

        logBox.gbc.fill = gbc.HORIZONTAL;
        logBox.addEndRow(logBean);
    	addEndRow(logBox);

        resampleBean = new StsBooleanFieldBean(wizard, "resampleLog", false, "Resample Log", false);
        intervalBean = new StsIntFieldBean(wizard, "numberOfSampleInterval", true, "Sample Interval:", true);
        computeMethodBean = new StsComboBoxFieldBean(this, "computeMethodString", "Compute Method:", computeMethodStrings);
        computeBox.gbc.fill = gbc.HORIZONTAL;
        computeBox.addToRow(resampleBean);
        computeBox.addEndRow(intervalBean);
        computeBox.addEndRow(computeMethodBean);
    	addEndRow(computeBox);

    	return true;
    }

    public void initialize()
    {
        // Build a list of unique logs from the selected wells
        Object[] wells = wizard.selectedWells;
        for(int i=0; i<wells.length; i++)
        {
            StsLogCurve[] curves = (StsLogCurve[])((StsWell)wells[i]).getLogCurves().getCastListCopy();
            // Determine if logs from this well are already in list.
            for(int j=0; j<curves.length; j++)
            {
                boolean inList = alreadyInList(curves[j]);
                if(!inList)
                    availableLogCurves = (StsLogCurve[])StsMath.arrayAddElement(availableLogCurves, curves[j]);
            }
        }
        logBean.initialize(this, "selectedLog", "Log:", availableLogCurves);

        // If valid depth grid already exists sample interval must be equal to that of grid.
        if(wizard.gridDefinitionExists())
        {
            resampleBean.setEditable(false);
            intervalBean.setEditable(false);
            intervalBean.setValue(1); // Sample in every cell.
            resampleBean.setSelected(true);  // Resample to project grid.
        }
        wizard.rebuild();
    }

    public boolean alreadyInList(StsLogCurve log)
    {
        if(availableLogCurves == null) return false;
        for(int k=0; k<availableLogCurves.length; k++)
        {
            if(availableLogCurves[k].getName().equalsIgnoreCase(log.getName()))
                return true;
        }
        return false;
    }
    
    public void setComputeMethodString(String method)
    {
        for(int i=0; i<computeMethodStrings.length; i++)
        {
            if(computeMethodStrings[i].equalsIgnoreCase(method))
                computeMethod = (byte)i;
        }
    }

    public String getComputeMethodString() { return computeMethodStrings[computeMethod]; }
    public byte getComputeMethod() { return computeMethod; }
    public StsLogCurve getSelectedLog() { return selectedLog; }
    public void setSelectedLog(StsLogCurve log) { selectedLog = log; }
}