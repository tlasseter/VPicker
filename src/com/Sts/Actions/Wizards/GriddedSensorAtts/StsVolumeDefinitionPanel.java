package com.Sts.Actions.Wizards.GriddedSensorAtts;

import com.Sts.DBTypes.*;
import com.Sts.MVC.StsProject;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.UI.Table.*;
import com.Sts.Actions.Wizards.PreStack3d.StsPreStackWizard;
import com.Sts.Actions.Wizards.PreStack3d.StsPreStackSurveyDefinition;
import com.Sts.Actions.Wizards.PreStack3d.Sts3PointSurveyDefinitionPanel;

import javax.swing.*;
import java.awt.*;
import java.util.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author T.Lasseter
 * @version 1.1
 */

public class StsVolumeDefinitionPanel extends StsJPanel
{
	private StsGriddedSensorAttsWizard wizard;
    private StsVolumeDefinition wizardStep;
    private StsSeismicBoundingBox boundingBox = new StsSeismicBoundingBox(false);

    private StsGroupBox surveyDefinitionBox;
	private StsJPanel parameterPanel;

	private StsDoubleFieldBean xMinBean;
	private StsDoubleFieldBean yMinBean;
	private StsFloatFieldBean zStartBean;
    private StsDoubleFieldBean xMaxBean;
	private StsDoubleFieldBean yMaxBean;
	private StsFloatFieldBean zEndBean;
	private StsFloatFieldBean inlineOriginBean;
	private StsFloatFieldBean xlineOriginBean;
	private StsFloatFieldBean xIncBean;
	private StsFloatFieldBean yIncBean;
	private StsFloatFieldBean zIncBean;
	private StsFloatFieldBean angleBean;

    public StsVolumeDefinitionPanel(StsGriddedSensorAttsWizard wizard, StsVolumeDefinition wizardStep)
	{
		super(true); // true adds insets
		this.wizard = wizard;
        this.wizardStep = wizardStep;
        try
		{
			constructBeans();
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
	}

	private void constructBeans()
	{
		xMinBean = new StsDoubleFieldBean(this, "xMin", true, "X Minimum:", true);
        xMinBean.setToolTipText("The minimum x coordinate inside the project bounds.");
        xMaxBean = new StsDoubleFieldBean(this, "xMax", true, "X Maximum:", true);
        xMinBean.setToolTipText("The maximum x coordinate inside the project bounds.");

		yMinBean = new StsDoubleFieldBean(this, "yMin", true, "Y Minimum:", true);
        xMinBean.setToolTipText("The minimum y coordinate inside the project bounds.");
        yMaxBean = new StsDoubleFieldBean(this, "yMax", true, "Y Maximum:", true);
        yMinBean.setToolTipText("The maximum y coordinate inside the project bounds.");

		zStartBean = new StsFloatFieldBean(boundingBox, "zMin", true, "Z Minimum:", true);
        zStartBean.setToolTipText("The minimum z coordinate inside the project bounds.");
        zEndBean = new StsFloatFieldBean(boundingBox, "zMax", true, "Z Maximum:", true);
        zEndBean.setToolTipText("The maximum z coordinate inside the project bounds.");

		inlineOriginBean = new StsFloatFieldBean(boundingBox, "rowNumMin", true, "InLine Origin:", true);
        inlineOriginBean.setToolTipText("The inline number corresponding to minimum Y");
		xlineOriginBean = new StsFloatFieldBean(boundingBox, "colNumMin", true, "Xline Origin:", true);
        xlineOriginBean.setToolTipText("The crossline number corresponding to minimum Y");

		xIncBean = new StsFloatFieldBean(boundingBox, "xInc", true, "X Interval:", true);
        xIncBean.setToolTipText("Grid increment in the X direction.");
		yIncBean = new StsFloatFieldBean(boundingBox, "yInc", true, "Y Interval:", true);
        yIncBean.setToolTipText("Grid increment in the Y direction.");
		zIncBean = new StsFloatFieldBean(boundingBox, "zInc", true, "Z Interval:", true);
        zIncBean.setToolTipText("Grid increment in the Z direction.");
		angleBean = new StsFloatFieldBean(boundingBox, "angle", 0.0f, 360.0f, "Angle:", true);
        angleBean.setValue(0.0f);
		angleBean.setEditable(false);
    }

	public void initialize()
	{
        float extent = wizard.getModel().getProject().getXMax() - wizard.getModel().getProject().getXMin();
        int nSteps = 100;
        int xyStepSize = (int)extent/nSteps;
        extent = wizard.getModel().getProject().getZorTMax() - wizard.getModel().getProject().getZorTMin();
        nSteps = 100;
        int zStepSize = (int)extent/nSteps;

        xMinBean.setValue(wizard.getModel().getProject().getXOrigin() + wizard.getModel().getProject().getXMin());
        yMinBean.setValue(wizard.getModel().getProject().getYOrigin() + wizard.getModel().getProject().getYMin());
        zStartBean.setValue(wizard.getModel().getProject().getZorTMin());
        xMaxBean.setValue(wizard.getModel().getProject().getXOrigin() + wizard.getModel().getProject().getXMax());
        yMaxBean.setValue(wizard.getModel().getProject().getYOrigin() + wizard.getModel().getProject().getYMax());
        zEndBean.setValue(wizard.getModel().getProject().getZorTMax());
        inlineOriginBean.setValue(1.0f);
        xlineOriginBean.setValue(1.0f);
        xIncBean.setValue(xyStepSize);
        yIncBean.setValue(xyStepSize);
        zIncBean.setValue(zStepSize);

        boundingBox.setXOrigin(wizard.getModel().getProject().getXOrigin());
        boundingBox.setYOrigin(wizard.getModel().getProject().getYOrigin());
        boundingBox.setXMin(wizard.getModel().getProject().getXMin());
        boundingBox.setYMin(wizard.getModel().getProject().getYMin());
        boundingBox.setZMin(wizard.getModel().getProject().getZorTMin());
        boundingBox.setXMax(wizard.getModel().getProject().getXMax());
        boundingBox.setYMax(wizard.getModel().getProject().getYMax());
        boundingBox.setZMax(wizard.getModel().getProject().getZorTMax());
        boundingBox.setXInc(xyStepSize);
        boundingBox.setYInc(xyStepSize);
        boundingBox.setZInc(zStepSize);
        boundingBox.setRowNumMin(1.0f);
        boundingBox.setColNumMin(1.0f);

        boundingBox.setZDomain(StsProject.TD_DEPTH);

        buildPanel();
	}

    private void buildPanel()
	{
        removeAll();

        surveyDefinitionBox = new StsGroupBox("Survey Definition");
	    parameterPanel = StsJPanel.addInsets();

        gbc.fill = GridBagConstraints.HORIZONTAL;

        parameterPanel.gbc.fill = GridBagConstraints.HORIZONTAL;
		parameterPanel.addEndRow(xMinBean);
        parameterPanel.addEndRow(xMaxBean);
        parameterPanel.addEndRow(xIncBean);

		parameterPanel.addEndRow(yMinBean);
		parameterPanel.addEndRow(yMaxBean);
        parameterPanel.addEndRow(yIncBean);

 		parameterPanel.addEndRow(zStartBean);
        parameterPanel.addEndRow(zEndBean);
		parameterPanel.addEndRow(zIncBean);

		parameterPanel.addEndRow(inlineOriginBean);
		parameterPanel.addEndRow(xlineOriginBean);
		parameterPanel.addEndRow(angleBean);

        surveyDefinitionBox.gbc.fill = GridBagConstraints.HORIZONTAL;
		surveyDefinitionBox.add(parameterPanel);

        this.add(surveyDefinitionBox);
        
        wizard.rebuild();
    }

    public void setXMin(double value)
    {
        boundingBox.setXMin(0.0f);
        boundingBox.setXOrigin(value);
    }
    public void setYMin(double value)
    {
        boundingBox.setYMin(0.0f);
        boundingBox.setYOrigin(value);
    }
    public void setXMax(double value)
    {
        boundingBox.setXMax((float)(value - boundingBox.getXOrigin()));
    }
    public void setYMax(double value)
    {
        boundingBox.setYMax((float)(value - boundingBox.getYOrigin()));
    }

    public double getXMin()
    {
        return (float)boundingBox.getXOrigin();
    }
    public double getYMin()
    {
        return (float)boundingBox.getYOrigin();
    }
    public double getXMax()
    {
        return boundingBox.getXMax() + (float)boundingBox.getXOrigin();
    }
    public double getYMax()
    {
        return boundingBox.getYMax() + (float)boundingBox.getYOrigin();
    }

    public StsSeismicBoundingBox getSeismicBoundingBox()
    {
        return boundingBox;
    }

    public void accept()
    {
       // new StsMessage(null, StsMessage.INFO, "All fields must be specified. Fill in missing fields or cancel.");
        return;
    }


    public void cancel()
    {
        return;
    }
}