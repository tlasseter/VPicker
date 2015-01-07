package com.Sts.Actions.Wizards.Monitor;

import com.Sts.Actions.Wizards.StsWizard;
import com.Sts.Actions.Wizards.StsWizardStep;
import com.Sts.UI.Beans.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsDefineWellPanel extends StsJPanel
{
    private StsMonitorWizard wizard;
    private StsDefineWell wizardStep;

    String name = "Well";
    double xOrigin = 0.0f, yOrigin = 0.0f, zOrigin = 0.0f;

    StsGroupBox nameBox = new StsGroupBox("");
    StsStringFieldBean nameBean = null;

    StsGroupBox defineBox = new StsGroupBox("Define Coordinate");
    StsDoubleFieldBean xOriginBean = null;
    StsDoubleFieldBean yOriginBean = null;
    StsDoubleFieldBean zOriginBean = null;

    public StsDefineWellPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsMonitorWizard)wizard;
        this.wizardStep = (StsDefineWell)wizardStep;

        try
        {
            constructBeans();
            jbInit();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void constructBeans()
    {
        nameBean = new StsStringFieldBean(this, "name", true, "Name:");

        xOriginBean = new StsDoubleFieldBean(this, "xOrigin", true, "X Origin:");
        yOriginBean = new StsDoubleFieldBean(this, "yOrigin", true, "Y Origin:");
        zOriginBean = new StsDoubleFieldBean(this, "zOrigin", true, "Z Origin:");
    }

    public void initialize()
    {
        xOrigin = wizard.getModel().getProject().getXOrigin();
        xOriginBean.getValueFromPanelObject();
        yOrigin = wizard.getModel().getProject().getYOrigin();
        yOriginBean.getValueFromPanelObject();
        zOrigin = wizard.getModel().getProject().getZorTMin();
        zOriginBean.getValueFromPanelObject();
    }

    void jbInit() throws Exception
    {
        gbc.fill = gbc.HORIZONTAL;
        gbc.anchor = gbc.NORTH;        
        nameBox.addEndRow(nameBean);

        defineBox.gbc.fill = gbc.HORIZONTAL;
        defineBox.gbc.anchor = gbc.NORTH;
        defineBox.gbc.gridwidth = 1;
        defineBox.addEndRow(xOriginBean);
        defineBox.addEndRow(yOriginBean);
        defineBox.addEndRow(zOriginBean);

        add(nameBox);
        add(defineBox);
    }

    public String getName() { return name; }
    public void setName(String name) { this.name = name; }

    public double getXOrigin() { return xOrigin; }
    public void setXOrigin(double x) { xOrigin = x; }
    public double getYOrigin() { return yOrigin; }
    public void setYOrigin(double y) { yOrigin = y; }
    public double getZOrigin() { return zOrigin; }
    public void setZOrigin(double z) { zOrigin = z; }

}