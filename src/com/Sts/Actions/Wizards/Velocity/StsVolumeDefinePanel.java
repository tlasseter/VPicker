package com.Sts.Actions.Wizards.Velocity;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import java.awt.*;

/**
 * <p>Title: jS2S development</p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author T.J.Lasseter
 * @version c51c
 */
public class StsVolumeDefinePanel extends StsJPanel
{
    StsGroupBox editPanel;
    StsComboBoxFieldBean volumeTypeBean = null;
    StsFloatFieldBean minVelocityBean;
    StsFloatFieldBean maxVelocityBean;
    StsComboBoxFieldBean velocityUnitsBean;
    StsComboBoxFieldBean oneTwoWayVelBean;
    StsSeismicVolume volume = null;
    StsVelocityWizard wizard = null;
    StsWizardStep wizardStep;
    String velocityUnits = StsParameters.VEL_UNITS_NONE;
    String oneOrTwoWayVelocity = StsParameters.ONE_WAY_VELOCITY;

    double scaleMultiplier = 1.0f;

    static final String noneString = new String("none");

    public StsVolumeDefinePanel(StsVolumeDefine volumeDefine)
    {
		super();
        wizard = (StsVelocityWizard) volumeDefine.getWizard();
        wizardStep = volumeDefine;
        constructPanel();
        defineVelocityUnits();
    }

    private void defineVelocityUnits()
    {
        StsProject project = wizard.getModel().getProject();
        velocityUnits = project.getVelocityUnits();
    }

    private void constructPanel()
    {
        editPanel = new StsGroupBox("Define Selected PostStack3d");
        volumeTypeBean = new StsComboBoxFieldBean(wizard, "velocityTypeString", "Velocity type:", StsVelocityWizard.VELOCITY_TYPE_STRINGS);
        minVelocityBean = new StsFloatFieldBean(wizard, "minVelocity", 0.0f, StsParameters.largeFloat, "Min velocity:");
        maxVelocityBean = new StsFloatFieldBean(wizard, "maxVelocity", 0.0f, StsParameters.largeFloat, "Max velocity:");
        velocityUnitsBean = new StsComboBoxFieldBean(this, "velocityUnits", "Velocity units:", StsParameters.VEL_UNITS );
        oneTwoWayVelBean = new StsComboBoxFieldBean(this, "oneOrTwoWayVelocity", "One- or two-way velocity:", StsParameters.ONE_OR_TWO_WAY );
        editPanel.add(volumeTypeBean);
        editPanel.add(minVelocityBean);
        editPanel.add(maxVelocityBean);
        editPanel.add(velocityUnitsBean);
        editPanel.add(oneTwoWayVelBean);
        volume = wizard.getVelocityVolume();
        if(volume != null)
             setVolumeValues();
        this.add(editPanel,  new GridBagConstraints(0, 0, 1, 1, 1.0, 0.0
                ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(3, 2, 0, 0), 0, 0));
    }

    public void initializePanel()
    {
        volume = wizard.getVelocityVolume();
        if(volume != null)
        {
            editPanel.setLabel("Define " + wizard.getVelocityVolume().getName());
            setVolumeValues();
        }
    }

    private void setVolumeValues()
    {
        if (volume == null)
        {
            minVelocityBean.setValue(0.0f);
            maxVelocityBean.setValue(0.0f);
            velocityUnitsBean.setSelectedItem(StsParameters.VEL_UNITS_NONE);
        }
        else
        {
            float dataMin = volume.getDataMin();
            minVelocityBean.setValue(dataMin);
            float dataMax = volume.getDataMax();
            maxVelocityBean.setValue(dataMax);
            wizard.setMinVelocity(dataMin);
            wizard.setMaxVelocity(dataMax);

            // guess the velocity units and whether one or two way
            oneOrTwoWayVelocity = StsParameters.ONE_WAY_VELOCITY;
            velocityUnits = StsParameters.estimateVelocityUnits(dataMin);
            if(dataMin < 1.0f)
                oneOrTwoWayVelocity = StsParameters.TWO_WAY_VELOCITY;

            velocityUnitsBean.setSelectedItem(velocityUnits);
            oneTwoWayVelBean.setSelectedItem(oneOrTwoWayVelocity);
        }
    }

    public void setVelocityUnits(String unitsString)
    {
        velocityUnits = unitsString;
        setVelScaleMultiplier();
    }

    public String getVelocityUnits() { return velocityUnits; }

    public void setOneOrTwoWayVelocity(String oneOrTwoWay)
    {
        this.oneOrTwoWayVelocity = oneOrTwoWay;
        setVelScaleMultiplier();
    }

    public String getOneOrTwoWayVelocity() { return oneOrTwoWayVelocity; }

    private void setVelScaleMultiplier()
    {
        StsProject project = wizard.getModel().getProject();
        scaleMultiplier = project.calculateVelScaleMultiplier(velocityUnits);

        if(oneOrTwoWayVelocity == StsParameters.ONE_WAY_VELOCITY)
            scaleMultiplier /= 2.0f;
    }

    public double getVelScaleMultiplier() { return scaleMultiplier; }
}
