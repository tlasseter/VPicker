package com.Sts.Actions.Wizards.Monitor;

import com.Sts.Actions.Wizards.*;
import com.Sts.UI.Beans.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsDefineSensorPanel extends StsJPanel
{
    private StsMonitorWizard wizard;
    private StsDefineSensor wizardStep; 
    
    String name = "Sensor";
    boolean isStatic = false;
    boolean isRelative = false;
    double xOrigin = 0.0f, yOrigin = 0.0f, zOrigin = 0.0f;
    boolean hasDate = true;
    String startDate = "";
    
    StsGroupBox nameBox = new StsGroupBox("");    
    StsStringFieldBean nameBean = null;    
    
    StsGroupBox defineBox = new StsGroupBox("Define Coordinate");
    StsBooleanFieldBean isStaticBean = null;  
    StsBooleanFieldBean isRelativeBean = null;    
    StsDoubleFieldBean xOriginBean = null;
    StsDoubleFieldBean yOriginBean = null;
    StsDoubleFieldBean zOriginBean = null;
    
    StsGroupBox timeDateBox = new StsGroupBox("Define Time-Date");    
    StsBooleanFieldBean hasDateBean = null;
    StsStringFieldBean startDateBean = null;
    
    public StsDefineSensorPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsMonitorWizard)wizard;
        this.wizardStep = (StsDefineSensor)wizardStep;

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
        
    	isStaticBean = new StsBooleanFieldBean(this, "isStatic", false, "Sensor Location is Static");
    	isRelativeBean = new StsBooleanFieldBean(this, "isRelative", false, "Sensor Locations are Relative");
        xOriginBean = new StsDoubleFieldBean(this, "xOrigin", true, "X Origin:");
        yOriginBean = new StsDoubleFieldBean(this, "yOrigin", true, "Y Origin:");
        zOriginBean = new StsDoubleFieldBean(this, "zOrigin", true, "Z Origin:");
    	
        hasDateBean = new StsBooleanFieldBean(this, "hasDate", true, "Input will have a Date Field");
        startDateBean = new StsStringFieldBean(this, "startDate", true, "Start Date:");        
    }

    public void initialize()
    {
        xOrigin = wizard.getModel().getProject().getXOrigin();
        xOriginBean.getValueFromPanelObject();
        yOrigin = wizard.getModel().getProject().getYOrigin();
        yOriginBean.getValueFromPanelObject();        
        zOrigin = wizard.getModel().getProject().getZorTMin();
        zOriginBean.getValueFromPanelObject();
		startDateBean.setEditable(false);
		hasDateBean.setValue(true);
		isRelativeBean.setValue(false);
        setIsRelative(false);
		isStaticBean.setValue(false);
		configureCoordinateComponents();
    }

    void jbInit() throws Exception
    {
        gbc.fill = gbc.HORIZONTAL;
        gbc.anchor = gbc.WEST;        
        nameBox.addEndRow(nameBean);

        defineBox.gbc.fill = gbc.HORIZONTAL;
        defineBox.gbc.anchor = gbc.WEST;
        defineBox.gbc.gridwidth = 2;
        defineBox.addEndRow(isStaticBean);
        defineBox.addEndRow(isRelativeBean);
        defineBox.gbc.gridwidth = 1;
        defineBox.addEndRow(xOriginBean);
        defineBox.addEndRow(yOriginBean);
        defineBox.addEndRow(zOriginBean);
        
        timeDateBox.gbc.fill = gbc.HORIZONTAL; 
        timeDateBox.gbc.anchor = gbc.WEST; 
        timeDateBox.gbc.gridwidth = 2;
        timeDateBox.addEndRow(hasDateBean);
        timeDateBox.gbc.gridwidth = 1;
        timeDateBox.addEndRow(startDateBean);
        
        add(nameBox);
        add(defineBox);
        add(timeDateBox);
    }

    public String getName() { return name; }
    public void setName(String name) { this.name = name; }
    public boolean getIsStatic() { return isStatic; }
    public void setIsStatic(boolean val) 
    { 
    	isStatic = val;
    	if(isStatic)
    	{
    		isRelativeBean.setValue(false);
    		isRelativeBean.setEditable(false);
    	}
    	else
    		isRelativeBean.setEditable(true);   	    	
    	configureCoordinateComponents();
    }
    public void configureCoordinateComponents()
    {
    	if(isStatic)
    	{
    		xOriginBean.setEditable(true);
    		yOriginBean.setEditable(true);
    		zOriginBean.setEditable(true);
    	}
    	else
    	{
    		if(isRelative)
    		{
        		xOriginBean.setEditable(true);
        		yOriginBean.setEditable(true);
        		zOriginBean.setEditable(true);    			
    		}
    		else
    		{
        		xOriginBean.setEditable(false);
        		yOriginBean.setEditable(false);
        		zOriginBean.setEditable(false);     			
    		}
    	}    	
    }
    public boolean getIsRelative() { return isRelative; }
    public void setIsRelative(boolean val) 
    { 
    	isRelative = val; 
    	if(isRelative)
    	{
    		isStaticBean.setValue(false);
    		isStaticBean.setEditable(false);
    	}
    	else
    		isStaticBean.setEditable(true);
    	configureCoordinateComponents();
    }
    public double getXOrigin() { return xOrigin; }
    public void setXOrigin(double x) { xOrigin = x; }
    public double getYOrigin() { return yOrigin; }
    public void setYOrigin(double y) { yOrigin = y; }
    public double getZOrigin() { return zOrigin; }
    public void setZOrigin(double z) { zOrigin = z; }
    
    public boolean getHasDate() { return hasDate; }
    public void setHasDate(boolean val) 
    { 
    	hasDate = val; 
    	if(!hasDate)
    		startDateBean.setEditable(true);
    	else
    		startDateBean.setEditable(false);
    }
    public String getStartDate() { return startDate; }    
    public void setStartDate(String date) 
    { 
    	startDate = date; 
    }
}
