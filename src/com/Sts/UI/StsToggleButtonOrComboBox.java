package com.Sts.UI;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.UI.Icons.*;
import com.Sts.UI.Toolbars.*;

import javax.swing.*;
import java.awt.*;
import java.util.*;

public class StsToggleButtonOrComboBox
{
    StsToolbar toolbar = null;
    String componentName = null;
    StsToggleButtonComboBox comboBox = new StsToggleButtonComboBox();
    StsToggleButton toggleButton = null;

    ArrayList buttons = new ArrayList(4);

    public StsToggleButtonOrComboBox()
    {
    }

    public StsToggleButtonOrComboBox(StsToolbar toolbar)
    {
        this.toolbar = toolbar;
    }

    public void addButton(StsToggleButton button)
    {
        buttons.add(button);
        comboBox.addButton(button);

        int nButtons = buttons.size();
        if(nButtons == 2)
            replaceComponent(toggleButton, comboBox);
        else if(nButtons == 1)
        {
            toggleButton = button;
            addComponent(toggleButton);
        }
    }

    public void deleteButton(String name)
    {
        comboBox.deleteButton(name);
        int nButtons = buttons.size();
        if(nButtons == 1)
            toggleButton = (StsToggleButton)buttons.get(0);
        else if(nButtons == 0)
            removeComponent(toggleButton);
    }

    private void addComponent(Component component)
    {
        if(toolbar == null) return;
        toolbar.addComponent(component);
    }

    private void replaceComponent(Component oldComponent, Component newComponent)
    {
        if(toolbar == null) return;
        toolbar.replaceComponent(oldComponent, newComponent);
    }

    private void removeComponent(Component component)
    {
        if(toolbar == null) return;
        toolbar.removeComponent(component);
    }

    public void setSelectedButton(StsToggleButton button, boolean isSelected)
    {
        comboBox.setSelectedButton(button, isSelected);
    }

    static public void main(String[] args)
    {
        try
        {
            JDialog dialog = new JDialog((Frame)null, "Button Panel Test", true);
            StsToolbar toolbar = new StsToolbar("Toolbar Test");
            dialog.getContentPane().add(toolbar);
            StsToggleButtonOrComboBox buttonBox = new StsToggleButtonOrComboBox(toolbar);
            StsToggleButtonOrComboBoxTest test = new StsToggleButtonOrComboBoxTest();

            StsToggleButton seismicButton = new StsToggleButton("Seismic", "Seismic tip.", test, "doSeismic");
            seismicButton.setText("Seismic");
            seismicButton.setIcon(StsIcon.createIcon("seismic.gif"));
            buttonBox.addButton(seismicButton);

            StsToggleButton surfacesButton = new StsToggleButton("Surfaces", "Surfaces tip.", test, "doSurfaces");
            surfacesButton.setText("Surfaces");
            surfacesButton.setIcon(StsIcon.createIcon("surfaces.gif"));
            buttonBox.addButton(surfacesButton);

            dialog.pack();
            dialog.setVisible(true);
        }
        catch(Exception e)
        {
            e.printStackTrace();
        }
    }
}

class StsToggleButtonOrComboBoxTest
{
    StsToggleButtonOrComboBoxTest()
    {
    }

    public void doSeismic()
    {
        System.out.println("Do seismic.");
    }

    public void doSurfaces()
    {
        System.out.println("Do surfaces.");
    }
}

