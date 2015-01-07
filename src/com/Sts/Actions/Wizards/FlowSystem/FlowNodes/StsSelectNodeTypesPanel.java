package com.Sts.Actions.Wizards.FlowSystem.FlowNodes;

import com.Sts.Actions.Wizards.*;
import com.Sts.UI.Beans.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSelectNodeTypesPanel extends JPanel implements ActionListener
{
    private StsWizard wizard;
    private StsWizardStep wizardStep;

    StsGroupBox typeBox = new StsGroupBox("Types (selected two)");

    JRadioButton tanksBtn = new JRadioButton("Tanks");
    JRadioButton pumpsBtn = new JRadioButton("Pumps");
    JRadioButton valvesBtn = new JRadioButton("Valves");
    JRadioButton pipelinesBtn = new JRadioButton("Pipelines");
    
    byte TANKS = 0;
    byte PUMPS = 1;
    byte VALVES = 2;
    byte PIPELINES = 3;

    byte type = TANKS;
    
    JTextArea typeDescription = new JTextArea();
    ButtonGroup typeGroup = new ButtonGroup();

    GridBagLayout gridBagLayout1 = new GridBagLayout();
    GridBagLayout gridBagLayout2 = new GridBagLayout();

    JPanel panel = new JPanel(new GridBagLayout());

    public StsSelectNodeTypesPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = wizardStep;
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
    	tanksBtn.setEnabled(true);
    	pumpsBtn.setEnabled(true);
    	pipelinesBtn.setEnabled(false);
    	valvesBtn.setEnabled(false);	
    }
    void jbInit() throws Exception
    {
        this.setLayout(gridBagLayout2);

        typeBox.gbc.anchor = typeBox.gbc.WEST;
        typeBox.add(tanksBtn);
        typeBox.add(pumpsBtn);
        typeBox.add(valvesBtn);
        typeBox.add(pipelinesBtn);

        tanksBtn.addActionListener(this);
        pumpsBtn.addActionListener(this);
        valvesBtn.addActionListener(this);
        pipelinesBtn.addActionListener(this);

        typeDescription.setBackground(Color.lightGray);
        typeDescription.setFont(new java.awt.Font("Dialog", 0, 10));
        typeDescription.setBorder(BorderFactory.createEtchedBorder());
        typeDescription.setText("Storage tanks which may be limited by valve, pump or pipeline capacities.");
        typeDescription.setLineWrap(true);
        typeDescription.setWrapStyleWord(true);
        panel.add(typeBox,  new GridBagConstraints(0, 0, 1, 2, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.VERTICAL, new Insets(0, 9, 0, 0), 0, 0));
        panel.add(typeDescription,  new GridBagConstraints(1, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 5, 0, 5), 0, 0));
        this.add(panel,  new GridBagConstraints(0, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0));
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();
        if(source == pumpsBtn)
        {
            typeDescription.setText("Pumps which may be limited by valve, tank or pipeline capacities.");
            type = PUMPS;
        }
        if(source == tanksBtn)
        {
            typeDescription.setText("Storage tanks which may be limited by valve, pump or pipeline capacities.");
            type = TANKS;
        }
        if(source == pipelinesBtn)
        {
            typeDescription.setText("Pipelines which may be limited by valve or pump capacities.");
            type = PIPELINES;
        }
        if(source == valvesBtn)
        {
            typeDescription.setText("Valves which may be limited tank, pipeline or pump capacities.");
            type = VALVES;
        }
    }

    public byte[] getNodeTypes()
    {
    	byte[] types = new byte[2];
    	int cnt = 0;
    	if(tanksBtn.isSelected())
    	{
    		types[cnt] = TANKS;
    		cnt++;
    	}
    	if(pipelinesBtn.isSelected())
    	{
    		types[cnt] = PIPELINES;
    		cnt++;
    	}
    	if(valvesBtn.isSelected())
    	{
    		types[cnt] = VALVES;
    		cnt++;
    	}
    	if(pumpsBtn.isSelected())
    	{
    		types[cnt] = PUMPS;
    		cnt++;
    	}
        return types;
    }

}
