package com.Sts.Actions.Wizards.AncillaryData;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2002</p>
 *
 * <p>Company: </p>
 *
 * @author not attributable
 * @version 1.0
 */

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;

import java.awt.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsAncillaryDataDefinePanel extends StsJPanel
{
    private StsAncillaryDataWizard wizard;
    private StsAncillaryDataDefine wizardStep;
    private StsModel model = null;
    private StsAncillaryData selectedObject = null;
    private StsColor stsColor;

    StsGroupBox box = new StsGroupBox();
    StsListFieldBean fileList = new StsListFieldBean();
    StsGroupBox parametersBox = new StsGroupBox("Selected File Parameters");
    StsButtonFieldBean updateBean = new StsButtonFieldBean();

    StsStringFieldBean commandStringBean = new StsStringFieldBean();
    StsStringFieldBean osNameBean = new StsStringFieldBean();
    StsDoubleFieldBean xLocBean = new StsDoubleFieldBean();
    StsDoubleFieldBean yLocBean = new StsDoubleFieldBean();
    StsDoubleFieldBean zLocBean = new StsDoubleFieldBean();
    StsCheckbox assignToWellCheckbox = new StsCheckbox("Assign to Well", "Do you want to tie any of the files to wells?");
    StsColorComboBoxFieldBean colorComboBoxBean = new StsColorComboBoxFieldBean();

    public StsAncillaryDataDefinePanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsAncillaryDataWizard)wizard;
        this.wizardStep = (StsAncillaryDataDefine)wizardStep;
        model = wizard.getModel();
        try
        {
            StsProject project = model.getProject();

            jbInit();
        }
        catch (Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initialize()
    {
        box.gbc.fill = gbc.HORIZONTAL;
        StsAncillaryData[] objects = wizard.getAncillaryDataObjects();
        if(objects.length > 0)
        {
            fileList.initialize(this, "selectedObject", "", objects);
        }
        box.addToRow(fileList);
        box.addEndRow(parametersBox);

        StsProject project = model.getProject();
        updateBean.initialize("Update", "Store new parameters", this, "fileParameters");
        StsColor[] colors = model.getSpectrum("Basic").getStsColors();
		stsColor = colors[0];
        colorComboBoxBean.initializeColors(this, "stsColor", "Color:", colors);
        commandStringBean.initialize(this, "commandString", true, "Launch Command:");
        osNameBean.initialize(this, "osAsString", false, "for Operating System:");
        String osName = System.getProperty("os.name").toString();
        osNameBean.setValue(osName);

        double[] xy0 = project.getAbsoluteXYCoordinates(project.getXMin(), project.getYMin());
        double[] xy1 = project.getAbsoluteXYCoordinates(project.getXMax(), project.getYMax());
        xLocBean.initialize(this, "xLoc", xy0[0], xy1[0], "X Location:");
        yLocBean.initialize(this, "yLoc", xy0[1], xy1[1], "Y Location:");
        zLocBean.initialize(this, "zLoc", project.getZorTMin(), project.getZorTMax(), "Z Location:");
        parametersBox.addEndRow(commandStringBean);
        parametersBox.addEndRow(osNameBean);
        parametersBox.addEndRow(colorComboBoxBean);
        parametersBox.addEndRow(xLocBean);
        parametersBox.addEndRow(yLocBean);
        parametersBox.addEndRow(zLocBean);
        parametersBox.addEndRow(updateBean);
        gbc.fill = gbc.HORIZONTAL;
        add(box);
        assignToWellCheckbox.setSelected(false);
        add(assignToWellCheckbox);
        validate();
    }

    public StsAncillaryData getSelectedObject()
    {
        return selectedObject;
    }
    public void setSelectedObject(StsAncillaryData object)
    {
        selectedObject = object;
        updateBeans();
    }

    public boolean getAssignToWell() { return assignToWellCheckbox.isSelected(); }
    public void updateBeans()
    {
        commandStringBean.setValue(selectedObject.getCommandString());
        osNameBean.setValue(selectedObject.getOsAsString());
        colorComboBoxBean.setValueObject(selectedObject.getStsColor());
        xLocBean.setValue(selectedObject.getXLoc());
        yLocBean.setValue(selectedObject.getYLoc());
        zLocBean.setValue(selectedObject.getZLoc());
    }

    public void fileParameters()
    {
        ;
    }

    void jbInit() throws Exception
    {
        this.setLayout(new GridBagLayout());
        add(box);
    }

    public double getXLoc()
    {
        if(selectedObject != null) return selectedObject.getXLoc();
        else return 0.0f;
    }

    public double getYLoc()
    {
        if(selectedObject != null) return selectedObject.getYLoc();
        else return 0.0f;
    }

    public double getZLoc()
    {
        if(selectedObject != null) return selectedObject.getZLoc();
        else return 0.0f;
    }

    public String getCommandString()
    {
        if(selectedObject != null) return selectedObject.getCommandString();
        else return null;
    }

    public void setXLoc(double value)
    {
        if(selectedObject != null)
            selectedObject.setXLoc(value);
    }

    public void setYLoc(double value)
    {
        if(selectedObject != null)
            selectedObject.setYLoc(value);
    }

    public void setZLoc(double value)
    {
        if(selectedObject != null)
            selectedObject.setZLoc(value);
    }

    public void setCommandString(String cmd)
    {
        if(selectedObject != null)
            selectedObject.setCommandString(cmd);
    }

    public void setOsAsString(String osName)
    {
        if(selectedObject != null)
            selectedObject.setOsAsString(osName);
    }
    public String getOsAsString()
    {
        if(selectedObject != null)
            return selectedObject.getOsAsString();
        return null;
    }

    public StsColor getStsColor()
    {
        if(selectedObject != null)
            return selectedObject.getStsColor();
        return null;
    }
    public void setStsColor(StsColor color)
    {
        if(selectedObject != null)
            selectedObject.setStsColor(color);
    }
}
