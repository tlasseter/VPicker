package com.Sts.Actions.Wizards.GriddedSensorAtts;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.StsSeismicBoundingBox;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.StsMath;
import com.Sts.Utilities.StsParameters;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSensorVVolumePanel extends StsJPanel
{
    private StsGriddedSensorAttsWizard wizard;
    private StsWizardStep wizardStep;

    private StsModel model = null;
    private String accumType = StsSensorVirtualVolume.ACCUM_METHOD[StsSensorVirtualVolume.ONOFF];
    Object selectedAttribute = null;
    byte shapeType = StsSensorVirtualVolume.SPHERE;
    float xyOffset = 100.0f;
    float zOffset = 100.0f;
    float zStep = 10.0f;
    float zSpan = 100.0f;

    StsGroupBox accumBox = new StsGroupBox("Define Accumulation Method");
    StsGroupBox searchBox = new StsGroupBox("Define Size & Shape");

    StsComboBoxFieldBean statMethodBean = new StsComboBoxFieldBean();
    StsComboBoxFieldBean attributeBean = new StsComboBoxFieldBean();
    StsComboBoxFieldBean shapeBean = new StsComboBoxFieldBean();
    private StsFloatFieldBean xyOffsetBean = new StsFloatFieldBean();
    private StsFloatFieldBean zOffsetBean = new StsFloatFieldBean();
    private StsFloatFieldBean zStepBean = new StsFloatFieldBean();
    private StsStringFieldBean nameBean = new StsStringFieldBean();

    public StsSensorVVolumePanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsGriddedSensorAttsWizard) wizard;
        this.wizardStep = wizardStep;
        this.model = wizard.getModel();
        try
        {
            jbInit();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initialize()
    {
        model = wizard.getModel();
        boolean hasCurves = false;

        Object[] sensors = (Object[])((StsGriddedSensorAttsWizard)wizard).getSelectedSensors();
        StsTimeCurve[] curves = ((StsSensor)sensors[0]).getPropertyCurvesExcludingXYZ();
        if(curves.length > 0)
            hasCurves = true;

        attributeBean.initialize(this,"attribute","Attribute:");
        StsTimeCurve curve = new StsTimeCurve(false);
        curve.setName("None");
        attributeBean.addItem(curve);
        for(int i=0; i<curves.length; i++)
            attributeBean.addItem(curves[i]);
        selectedAttribute = curve;

        shapeBean.initialize(this,"shape","Shape:",StsSensorVirtualVolume.shapes);

        float xExtent = 0;
        float[] zRange = new float[] {model.getProject().getZorTMin(), model.getProject().getZorTMin()};
        if(hasCurves)
            zRange = ((StsSensor)sensors[0]).getZTRange();
        float stepSize = (zRange[1] - zRange[0])/100.0f;
        zOffset = stepSize;
        zSpan = zRange[1] - zRange[0];
        if(((StsGriddedSensorAttsWizard) wizard).getBoundingBox() == null)
        {
            xExtent = model.getProject().getXMax() - model.getProject().getXMin();
            stepSize = xExtent/100;
            xyOffset = stepSize;
            zStepBean.setEditable(false);
        }
        else
        {
            StsSeismicBoundingBox box = ((StsGriddedSensorAttsWizard) wizard).getBoundingBox();
            xExtent = box.getXMax() - box.getXMin();
            stepSize = box.getXInc();
            xyOffset = box.getXInc();
            if(box.getZDomain() == StsParameters.TD_DEPTH)
                zStep = box.getZInc();
        }
        xyOffsetBean.setValueAndRangeFixStep(stepSize, xyOffset, xExtent, stepSize);
        zOffsetBean.setValueAndRangeFixStep(zStep, zStep, zSpan, zStep);
        zStepBean.setValueAndRangeFixStep(zStep, zStep, zSpan, zStep);
        zOffsetBean.setEditable(false);
    }

    void jbInit() throws Exception
    {
    	statMethodBean.initialize(this, "accumTypeString", "Accumulation Type:", StsSensorVirtualVolume.ACCUM_METHOD);
        statMethodBean.setToolTipText("Select the method used to compute the cell values.");
        zStepBean.initialize(this,"zStep",true,"Z Interval:",true);
        zStepBean.setToolTipText("Specify the Z interval for the cube to be computed.");
    	//attributeBean.initialize(this, "attribute", "Attribute:", new Object[0]);
        attributeBean.setToolTipText("Select the attribute to use.");
        xyOffsetBean.initialize(this,"xyOffset",true,"Radius:",true);
        xyOffsetBean.setToolTipText("Specify the horizontal radius to collect events for computation.");
        zOffsetBean.initialize(this,"zOffset",true,"Time/Height:",true);
        zOffsetBean.setToolTipText("Specify the vertical time/distance to collect events for computation.");
        nameBean.initialize(wizard, "volumeName", true, ((StsGriddedSensorAttsWizard)wizard).getVolumeName());
        nameBean.setToolTipText("Specify the name of the new virtual volume.");

        accumBox.addEndRow(nameBean);        
    	accumBox.addEndRow(statMethodBean);
    	accumBox.addEndRow(attributeBean);
        accumBox.addEndRow(zStepBean);

        searchBox.addEndRow(shapeBean);
        searchBox.addEndRow(xyOffsetBean);
        searchBox.addEndRow(zOffsetBean);

    	gbc.fill = gbc.HORIZONTAL;
    	add(accumBox);
        add(searchBox);
    }

    public void setAccumTypeString(String aType)
    {
        accumType = aType;       
    }
    public String getAccumTypeString() { return accumType; }
    public byte getAccumType()
    {
        return StsSensorVirtualVolume.getAccumMethodFromString(accumType);
    }
    public float getXyOffset()
    {
        return xyOffset;
    }
    public void setZOffset(float off)
    {
        zOffset = off;
    }
    public float getZOffset()
    {
        return zOffset;
    }
    public void setShape(String shape)
    {
        for(int i=0; i<StsSensorVirtualVolume.shapes.length; i++)
        {
            if(StsSensorVirtualVolume.shapes[i].equalsIgnoreCase(shape))
                shapeType = (byte)i;
        }
        if(shapeType == StsSensorVirtualVolume.SPHERE)
            zOffsetBean.setEditable(false);
        else
            zOffsetBean.setEditable(true);
    }
    public String getShape()
    {
        return StsSensorVirtualVolume.shapes[shapeType];
    }
    public byte getShapeType()
    {
        return shapeType;
    }
    public void setXyOffset(float off)
    {
        xyOffset = off;
    }
    public void setZStep(float step)
    {
        zStep = step;
        zOffsetBean.setValueAndRangeFixStep(zStep, zStep, zSpan, zStep);
    }
    public float getZStep() { return zStep; }
    public void setAttribute(Object attribute)
    {
        selectedAttribute = attribute;       
    }
    public Object getAttribute()
    {
        return selectedAttribute; 
    }
}
