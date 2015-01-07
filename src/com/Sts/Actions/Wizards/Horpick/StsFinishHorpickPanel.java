package com.Sts.Actions.Wizards.Horpick;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;

import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version beta 1.0
 */

public class StsFinishHorpickPanel extends StsFieldBeanPanel
{
    StsHorpickWizard wizard;
    StsGroupBox displayBox = new StsGroupBox("Display controls");
    StsComboBoxFieldBean propertyList = new StsComboBoxFieldBean();
    StsColorscalePanel corCoefsColorscalePanel = null;
    StsSliderBean minCorFilterSlider = new StsSliderBean();

    public StsFinishHorpickPanel(StsWizard wizard)
    {
        this.wizard = (StsHorpickWizard)wizard;

        try
        {
            constructPanel();
            propertyList.setValueObject(this.wizard.getDisplayPropertyName());
        }
        catch(Exception e)
        {
            e.printStackTrace();
        }
    }

    private void constructPanel() throws Exception
    {
        setLayout(new GridBagLayout());

        displayBox.setLayout(new GridBagLayout());
        add(displayBox, new GridBagConstraints(0, 0, 1, 1, 1.0, 0.0,
            GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(10, 5, 10, 10), 0, 0));

        initializeMinCorSlider();
//        minCorSlider.setLabelFormats("0.000", "0.000");
        displayBox.add(minCorFilterSlider, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0,
            GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(10, 10, 10, 10), 0, 0));

        String[] propertyNames = wizard.getDisplayPropertyNames();
        propertyList.initialize(wizard, "displayPropertyName", "Display Property", propertyNames);
        displayBox.add(propertyList, new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0,
            GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(10, 10, 10, 5), 0, 0));
    }

    public void initialize()
    {
        corCoefsColorscaleEnabled(false);
        initializePropertyList();
        initializeCorCoefsColorscalePanel();
        initializeFilterSlider();
    }

    private void initializePropertyList()
    {
        wizard.initializeDisplayPropertyName();
        propertyList.setSelectedItem(StsHorpick.displayPropertyPatchColor);
    }

    private void initializeFilterSlider()
    {
        StsHorpick horpick = wizard.getHorpick();
        if(horpick == null) return;
        float minCorrel = horpick.getMinCorrel();
        float minCorFilter = wizard.getMinCorFilter();
        minCorFilterSlider.initSliderValues(minCorrel, 1.0f, 0.01f, minCorFilter);
    }

    private void initializeCorCoefsColorscalePanel()
    {
        StsHorpick horpick = wizard.getHorpick();
        if(horpick == null) return;
//        if(corCoefsColorscalePanel != null) return;
        corCoefsColorscalePanel = horpick.getCorCoefsColorscalePanel();
    }

    public void corCoefsColorscaleEnabled(boolean enabled)
    {
        if(enabled == false)
        {
            if(corCoefsColorscalePanel != null) displayBox.remove(corCoefsColorscalePanel);
        }
        else
        {
            initializeCorCoefsColorscalePanel();
            displayBox.add(corCoefsColorscalePanel, new GridBagConstraints(0, 2, 1, 1, 1.0, 1.0,
            GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0));
        }
        validate();
    }

    private void initializeMinCorSlider()
    {
        minCorFilterSlider.initSliderValues(0.0f, 1.0f, 0.01f, 0.0f);
        minCorFilterSlider.addChangeListener
        (
            new ChangeListener()
            {
                public void stateChanged(ChangeEvent e)
                {
                    minCorSliderChanged(e);
                }
            }
        );
        minCorFilterSlider.addItemListener
        (
            new ItemListener()
            {
                public void itemStateChanged(ItemEvent e)
                {
                    minCorSliderItemChanged(e);
                }
            }
        );
    }

    private void minCorSliderChanged(ChangeEvent e)
    {
        float value = minCorFilterSlider.getValue();
        wizard.minCorSliderChanged(value);
    }

    private void minCorSliderItemChanged(ItemEvent e)
    {
//        wizard.minCorSliderItemChanged(value);
    }
}
