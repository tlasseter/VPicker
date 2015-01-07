package com.Sts.Actions.Wizards.VolumeStimulated;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSelectFractureSetsPanel extends StsJPanel
{
    private StsVolumeStimulatedWizard wizard;
    private StsSelectFractureSets wizardStep;

    private Object[] availableFractureSets = null;
    private Object[] selectedFractureSets = null;
    private StsModel model = null;

    StsJPanel beanPanel = new StsJPanel();
    StsGroupBox listBox = new StsGroupBox("Available Fracture Sets");
    StsListFieldBean fractureListBean;
    
    StsGroupBox panelBox = new StsGroupBox("Define Stimulated Volume");
    StsBooleanFieldBean enableBean = null;
    StsFloatFieldBean radiusBean = new StsFloatFieldBean();
    StsFloatFieldBean cellSizeBean = new StsFloatFieldBean();
    
	StsGroupBox analyzeBox = new StsGroupBox("Analyze Volume");
	StsComboBoxFieldBean unitsBean = new StsComboBoxFieldBean();
	StsStringFieldBean msgBean = new StsStringFieldBean(false);
	StsButton runButton;
	
    StsProgressPanel progressPanel;
    
    public StsSelectFractureSetsPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
    	this.wizard = (StsVolumeStimulatedWizard)wizard;
    	this.wizardStep = (StsSelectFractureSets)wizardStep;
    	model  = wizard.getModel();
        buildPanel();
    }

    public void buildPanel()
    {    	
    	progressPanel = StsProgressPanel.constructorWithCancelButton();
    	     	
    	availableFractureSets = (Object[])model.getCastObjectList(StsFractureSet.class);    	
    	
    	fractureListBean = new StsListFieldBean(this, "fractureSets", null, availableFractureSets);
    	gbc.fill = gbc.BOTH;
    	listBox.gbc.fill = gbc.fill;
    	listBox.add(fractureListBean);
    	add(listBox);
    	
    	enableBean = new StsBooleanFieldBean(this, "isVisible", "Enable");
    	enableBean.setToolTipText("To include a fracture set in the volume computation it must be enabled.");
    	
    	radiusBean.initialize(this, "radius", 0.f, 10000.f, "Stimulation Radius:", true);
    	radiusBean.fixStep(10.0f);
    	radiusBean.setToolTipText("Specify the stimulated radius for the selected fracture sets.");
        
    	cellSizeBean.initialize(this, "cellSize", 0.0f, 1000.0f, "Cell Size:", true);
    	cellSizeBean.fixStep(1.0f);
    	cellSizeBean.setToolTipText("Specify the cell size within the stimulated volume.");
                
    	panelBox.gbc.fill = gbc.HORIZONTAL;
    	panelBox.addEndRow(enableBean);
    	panelBox.addEndRow(radiusBean);
    	panelBox.addEndRow(cellSizeBean);
    	gbc.gridwidth = 2;
        add(panelBox);

        analyzeBox.gbc.fill = gbc.HORIZONTAL;
        msgBean.setText("Push button to start analysis");
        StsJPanel msgPanel = new StsJPanel();
        msgPanel.gbc.fill = gbc.HORIZONTAL;
        msgPanel.addEndRow(msgBean);
        
        StsJPanel btnPanel = new StsJPanel();
        unitsBean.initialize(wizard, "units", "Units:");
        unitsBean.setListItems(wizard.unitStrings, wizard.unitStrings[0]);
        runButton = new StsButton("Run Analysis", "Compute the volume within stimulated area.", this, "analyzeFractureSets", progressPanel);
        btnPanel.gbc.fill = gbc.NONE;
        btnPanel.addToRow(runButton);
        btnPanel.gbc.fill = gbc.HORIZONTAL;
        btnPanel.addEndRow(unitsBean);
        
        analyzeBox.addEndRow(msgPanel);
        analyzeBox.addEndRow(progressPanel);        
        analyzeBox.addEndRow(btnPanel);
        add(analyzeBox);               
    }

    public void initialize()
    {
    	((StsFractureSetClass)model.getStsClass(StsFractureSet.class)).setDisplayTypeVolume();
    }

    public void analyzeFractureSets(StsProgressPanel panel)
    {
        progressPanel.resetProgressBar();
        progressPanel.setValue(0);
        wizard.analyzeFractureSets(panel);
    }

    public void setFractureSets(Object fs)
    {
    	selectedFractureSets = (Object[])fractureListBean.getSelectedObjects();
    	radiusBean.setValue(((StsFractureSet)selectedFractureSets[0]).getStimulatedRadius());
    	cellSizeBean.setValue(((StsFractureSetClass)model.getStsClass(StsFractureSet.class)).getGridSize());    	
    	enableBean.setSelected(((StsFractureSet)selectedFractureSets[0]).getIsVisible());
    }

    public Object getFractureSets()
    {
        return availableFractureSets[0];
    }
    
    public Object[] getSelectedFractureSets() { return selectedFractureSets; } 
    
    public void setRadius(float radius)
    {
    	if(selectedFractureSets == null) return;
    	if(selectedFractureSets.length < 1) return;
    	for(int i=0; i<selectedFractureSets.length; i++)
    		((StsFractureSet)selectedFractureSets[i]).setStimulatedRadius(radius);
    	model.viewObjectRepaint(this, selectedFractureSets[0]);
    }
    public float getRadius()
    {
    	if(selectedFractureSets == null) 
    		return ((StsFractureSetClass)model.getStsClass(StsFractureSet.class)).getDefaultStimulatedRadius();    	
    	if(selectedFractureSets.length < 1) 
    		return ((StsFractureSetClass)model.getStsClass(StsFractureSet.class)).getDefaultStimulatedRadius();
    	return ((StsFractureSet)selectedFractureSets[0]).getStimulatedRadius();    	
    }
    public void setCellSize(float size)
    {
    	((StsFractureSetClass)model.getStsClass(StsFractureSet.class)).setGridSize(size);
    	if(selectedFractureSets != null)
    		model.viewObjectRepaint(this, selectedFractureSets[0]);
    }
    public float getCellSize()
    { 
    	return ((StsFractureSetClass)model.getStsClass(StsFractureSet.class)).getGridSize();    	
    }
    public boolean getIsVisible()
    {
    	if(selectedFractureSets == null) return false;
    	if(selectedFractureSets.length < 1) return false;    	
    	return ((StsFractureSet)selectedFractureSets[0]).getIsVisible(); 
    }
    public void setIsVisible(boolean visible)
    {
    	if(selectedFractureSets == null) return;
    	if(selectedFractureSets.length < 1) return;
    	for(int i=0; i<selectedFractureSets.length; i++)
    		((StsFractureSet)selectedFractureSets[i]).setIsVisible(visible);	
    }
    public void setMessage(String msg)
    {
    	msgBean.setText(msg);
    }    
}
