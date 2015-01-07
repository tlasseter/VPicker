package com.Sts.Actions.Wizards.SurfaceCurvature;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsProcessSurfacePanel extends StsJPanel
{
    StsSurface surface;
    byte attributeType = StsSurfaceCurvatureAttribute.CURVPos;
    public int filterSize = 3;

    private StsSurfaceCurvatureWizard wizard;
    private StsProcessSurface wizardStep;
    private StsSurfaceCurvatureAttribute attribute;

    StsGroupBox attributeBox = new StsGroupBox("Select Curvature Attribute");
    StsGroupBox filterBox = new StsGroupBox("Select Analysis Window");
	StsGroupBox analyzeBox = new StsGroupBox("Analyze Curvature");
	StsEditableColorscaleFieldBean colorscaleBean = null;
    StsIntFieldBean filterSizeBean;
    StsComboBoxFieldBean curvAttributeBean;
	StsStringFieldBean msgBean = new StsStringFieldBean(false);
	StsButton runButton;
	StsButton saveButton;
    //StsButton exportButton;
    public StsProgressPanel progressPanel = StsProgressPanel.constructor(true, 5, 50);

    String attrString = StsSurfaceCurvatureAttribute.CURVPosString;
    
	
    public StsProcessSurfacePanel(StsSurfaceCurvatureWizard wizard, StsProcessSurface wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = wizardStep;
        buildPanel();
    }

    public void buildPanel()
    {
    	curvAttributeBean = new StsComboBoxFieldBean(this, "attributeString", "Attribute:", StsSurfaceCurvatureAttribute.CURV_ATTRIBUTE_STRINGS);
    	curvAttributeBean.setToolTipText("Select Curvature Attribute.");

        filterSizeBean = new StsIntFieldBean(this, "filterSize", true, "Filter Size:", true);
        filterSizeBean.setValueAndRangeFixStep(3, 3, 21, 2);
        filterSizeBean.setToolTipText("Specify analysis window size.");
        
        attributeBox.gbc.fill = GridBagConstraints.HORIZONTAL;
        attributeBox.add(curvAttributeBean);
        
        filterBox.gbc.fill = GridBagConstraints.HORIZONTAL;
        filterBox.add(filterSizeBean);

        analyzeBox.gbc.fill = GridBagConstraints.HORIZONTAL;
        
        StsJPanel btnPanel = StsJPanel.addInsets();
        runButton = new StsButton("Run Analysis", "Push button to start analysis.", this, "analyze");
        saveButton = new StsButton("Save to Model", "Save the attribute surface to the model.", this, "saveToModel");
        //exportButton = new StsButton("Export View", "Export the surface with computed attributes.", wizard, "exportView");
        saveButton.setEnabled(false);

        btnPanel.gbc.fill = GridBagConstraints.NONE;
        btnPanel.addToRow(runButton);
        btnPanel.addToRow(saveButton);
        //btnPanel.addToRow(exportButton);
        
        //analyzeBox.addEndRow(msgPanel);
        analyzeBox.addEndRow(btnPanel);

        StsJPanel mainPane = new StsJPanel();
        StsJPanel leftPane = new StsJPanel();
        JPanel rightPane = new StsJPanel();
        rightPane.setLayout(new BorderLayout());
        mainPane.setMinimumSize(new Dimension(400, 284));
        leftPane.gbc.fill = GridBagConstraints.HORIZONTAL;
        leftPane.gbc.anchor = GridBagConstraints.WEST;
        leftPane.addEndRow(attributeBox);
        leftPane.addEndRow(filterBox);
        leftPane.addEndRow(analyzeBox);

        colorscaleBean = new StsEditableColorscaleFieldBean(StsSurface.class, "colorscale");
        rightPane.add(colorscaleBean, BorderLayout.CENTER);
        rightPane.setPreferredSize(new Dimension(154, 284));

        mainPane.gbc.fill = GridBagConstraints.HORIZONTAL;
        mainPane.addToRow(leftPane);
        mainPane.addToRow(rightPane);
        
        gbc.fill = GridBagConstraints.BOTH;
        addEndRow(mainPane);
        addEndRow(progressPanel);
    }

    public void initialize()
    {
        surface = wizard.getSelectedSurface();
        initializeAttribute();
        initializeColorscale();
    }

    private void initializeAttribute()
    {
        if(surface.curvatureAttribute != null)
        {
            setCurvatureAttribute(surface.curvatureAttribute);
            return;
        }
        StsSurfaceCurvatureAttribute[] attributes = surface.getCurvatureAttributes();
        if(attributes.length == 0) return;
        setCurvatureAttribute(attributes[0]);
    }

    public void initializeColorscale()
    {
        if(surface == null) return;
        StsSurfaceCurvatureAttribute attribute = surface.curvatureAttribute;
        if (attribute != null)
        {
            colorscaleBean.setValueObject(attribute.getColorscale());
            colorscaleBean.setHistogram(attribute.getHistogram());
            colorscaleBean.setVisible(true);
        }
        else
            colorscaleBean.setVisible(false);
        validate();
    }
    
    public void updateColorscale()
    {
        StsToolkit.runLaterOnEventThread
        (
            new Runnable()
            {
                public void run()
                {
                    initializeColorscale();
                }
            }
        );
    }

    public String getAttributeString()
    {
    	return StsSurfaceCurvatureAttribute.CURV_ATTRIBUTE_STRINGS[attributeType];
    }

    public void setAttributeString(String string)
    {
        attributeType = StsParameters.getStringMatchByteIndex(StsSurfaceCurvatureAttribute.CURV_ATTRIBUTE_STRINGS, string);
        checkSetDisplay();
    }

    private void checkSetDisplay()
    {
        StsSurfaceCurvatureAttribute[] curvatureAttributes = surface.getCurvatureAttributes();
        if(curvatureAttributes == null) return;
        for(StsSurfaceCurvatureAttribute curvatureAttribute : curvatureAttributes)
        {
            if(curvatureAttribute.curveType != attributeType) continue;
            if(curvatureAttribute.filterSize != filterSize) continue;
            setCurvatureAttribute(curvatureAttribute);
            return;
        }
        colorscaleBean.setVisible(false);
        runButton.setEnabled(true);
        saveButton.setEnabled(false);
    }

    private void setCurvatureAttribute(StsSurfaceCurvatureAttribute curvatureAttribute)
    {
        surface.curvatureAttribute = curvatureAttribute;
        curvAttributeBean.doSetValueObject(curvatureAttribute.getAttributeString());
        filterSizeBean.setValue(curvatureAttribute.filterSize);
        updateColorscale();
        runButton.setEnabled(false);
        saveButton.setEnabled(!curvatureAttribute.isPersistent());
        surface.setNewSurfaceTexture(curvatureAttribute.surfaceTexture);
    }

    public boolean analyze()
    {
    	String attributeMessage = "Computing " + StsSurfaceCurvatureAttribute.CURV_ATTRIBUTE_STRINGS[attributeType];
        progressPanel.appendLine(attributeMessage + " on surface (" + wizard.getSelectedSurface().getName() + ") with a " + filterSize + "x" + filterSize + " filter.");
        Main.logUsage();
        
        runCreateAttribute();
        saveButton.setEnabled(true);
        runButton.setEnabled(false);
        return true;
    }
    
    public boolean saveToModel()
    {  
    	String attributeMessage = "Saving " + StsSurfaceCurvatureAttribute.CURV_ATTRIBUTE_STRINGS[attributeType];
    	progressPanel.appendLine(attributeMessage + " on surface (" + wizard.getSelectedSurface().getName() + ") to model.");
        StsSurface surface = wizard.getSelectedSurface();
        surface.saveCurvatureAttribute();
        saveButton.setEnabled(false);
        return true;
    }
    
    public void runCreateAttribute()
    {
        Runnable runCurvSurface = new Runnable()
        {
            public void run()
            {
            	StsSurface surface = wizard.getSelectedSurface();
                surface.createCurvatureAttribute(attributeType, filterSize, progressPanel);
                updateColorscale();
            }
        };

        StsToolkit.runRunnable(runCurvSurface);
    }

    public void setFilterSize(int size)
    {
        filterSize = size;
        checkSetDisplay();
    }
    public int getFilterSize() { return filterSize; }

}