package com.Sts.Actions.Wizards.Workflow;

import com.Sts.Actions.Wizards.*;
import com.Sts.MVC.*;
import com.Sts.Utilities.*;
import com.Sts.WorkflowPlugIn.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.lang.reflect.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsDefineWorkflowPanel extends JPanel implements ActionListener
{
    private StsWizard wizard;
    private StsWizardStep wizardStep;

    private StsModel model = null;
    JPanel jPanel1 = new JPanel();
    JComboBox typesCombo = new JComboBox();
    JLabel jLabel1 = new JLabel();
    JTextArea typeDescription = new JTextArea();
    GridBagLayout gridBagLayout2 = new GridBagLayout();
    GridBagLayout gridBagLayout1 = new GridBagLayout();

    Class[] plugInClasses = new Class[0];
    public StsWorkflowPlugIn selectedPlugIn = null;
    public boolean selectedPlugInStatus = true;

    public StsDefineWorkflowPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = wizardStep;
        try
        {
            initialize();
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
        setWorkflowSelectComboBoxItems();
    }
    private void setWorkflowSelectComboBoxItems()
    {
        try
        {
            String[] workflowPlugInNames = model.workflowPlugInNames;
            if(workflowPlugInNames == null) return;
            int nPlugIns = workflowPlugInNames.length;

            PlugInClass[] classes = new PlugInClass[nPlugIns];
            for(int n = 0; n < nPlugIns; n++)
                classes[n] = new PlugInClass(workflowPlugInNames[n]);

            for(int n = 0; n < nPlugIns; n++)
            {
                Class plugInClass = classes[n].getPlugInClass();
                typesCombo.addItem(classes[n].className);
                plugInClasses = (Class[])StsMath.arrayAddElement(plugInClasses, plugInClass);
            }

            int index = model.getWorkflowPlugInIndex();
            if (index < 0)
                return;
            typesCombo.setSelectedIndex(index);
            plugInSelected(index);
        }
        catch(Exception e)
        {
            StsException.outputException("StsWorkflowSelectPanel.setWorkflowSelectComboBoxItems() failed.",
                e, StsException.WARNING);
        }
    }

    void jbInit() throws Exception
    {
        jPanel1.setBorder(BorderFactory.createEtchedBorder());
        jPanel1.setLayout(gridBagLayout1);
        this.setLayout(gridBagLayout2);
        jLabel1.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel1.setText("Types:");
        typesCombo.addActionListener(this);
        typeDescription.setBackground(Color.lightGray);
        typeDescription.setFont(new java.awt.Font("SansSerif", 0, 10));
        typeDescription.setBorder(BorderFactory.createEtchedBorder());
        typeDescription.setLineWrap(true);
        typeDescription.setWrapStyleWord(true);

        jPanel1.add(jLabel1,   new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(9, 7, 0, 0), 15, 3));
        jPanel1.add(typesCombo,   new GridBagConstraints(1, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 6), 299, 4));
        jPanel1.add(typeDescription,   new GridBagConstraints(1, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 0, 6), 0, 29));
        this.add(jPanel1,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(6, 5, 5, 3), 3, 7));

    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();
        if(source == typesCombo)
        {
            int index = typesCombo.getSelectedIndex();
            plugInSelected(index);
        }
    }

    public void plugInSelected(int index)
    {
        try
        {
            Class plugInClass = plugInClasses[index];
            if(plugInClass.equals(selectedPlugIn)) return;
            Constructor constructor = plugInClass.getConstructor(new Class[0]);
            selectedPlugIn = (StsWorkflowPlugIn)constructor.newInstance(new Class[0]);
            selectedPlugInStatus = Main.workflowStatus[index];
            String description = selectedPlugIn.getDescription();
            typeDescription.setText(description);
            model.setWorkflowPlugIn(selectedPlugIn);
            model.win3d.rebuildWorkflow(selectedPlugIn);
        }
        catch(Exception e)
        {
            StsException.outputException("StsWorkflowSelectPanel.plugInSelected() failed.",
                                         e, StsException.WARNING);
        }
    }

    class PlugInClass
    {
        String plugInName;
        String className;
        String description;

        public PlugInClass(String plugInName)
        {
            try
            {
                this.plugInName = plugInName;
                Class plugInClass = getPlugInClass();
                Constructor constructor = plugInClass.getConstructor(new Class[0]);
                StsWorkflowPlugIn workflowPlugIn = (StsWorkflowPlugIn)constructor.newInstance(new Class[0]);
                description = workflowPlugIn.getDescription();
                className = plugInName.substring(plugInName.lastIndexOf(".")+1, plugInName.length());
            }
            catch(Exception e)
            {
                StsException.outputException("StsWorkflowSelectPanel.PlugInClass.constructor() failed.",
                    e, StsException.WARNING);
            }
        }

        public Class getPlugInClass()
        {
            try
            {
                ClassLoader classLoader = getClass().getClassLoader();
                return classLoader.loadClass(plugInName);
            }
            catch(Exception e)
            {
                StsException.outputException("StsWorkflowSelectPanel.getPLugInClass() failed for class: " + plugInName,
                    e, StsException.WARNING);
                return null;
            }

        }
    }

    public byte getType()
    {
        return (byte)typesCombo.getSelectedIndex();
    }

}
