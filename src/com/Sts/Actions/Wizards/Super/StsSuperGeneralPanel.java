package com.Sts.Actions.Wizards.Super;

import com.Sts.Actions.Wizards.*;

import javax.swing.*;
import javax.swing.border.*;
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

public class StsSuperGeneralPanel extends JPanel implements ActionListener
{
    private StsSuperWizard wizard;
    private StsSuperGeneral wizardStep;

    Border border1;
    JLabel jLabel1 = new JLabel();
    JTextArea jTextArea1 = new JTextArea();
    JRadioButton yesBtn = new JRadioButton();
    JRadioButton noBtn = new JRadioButton();
    ButtonGroup yesNoGroup = new ButtonGroup();
    GridBagLayout gridBagLayout1 = new GridBagLayout();

    public StsSuperGeneralPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsSuperWizard)wizard;
        this.wizardStep = (StsSuperGeneral)wizardStep;
        try
        {
            jbInit();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    void jbInit() throws Exception
    {
        border1 = BorderFactory.createCompoundBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED,Color.white,Color.white,new Color(178, 178, 178),new Color(124, 124, 124)),BorderFactory.createEmptyBorder(10,10,10,10));
        this.setLayout(gridBagLayout1);
        jLabel1.setFont(new java.awt.Font("Dialog", 1, 11));
        jLabel1.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel1.setText(wizardStep.getQuestion());
        jTextArea1.setBackground(SystemColor.menu);
        jTextArea1.setFont(new java.awt.Font("SansSerif", 0, 11));
        jTextArea1.setBorder(BorderFactory.createEtchedBorder());
        jTextArea1.setMinimumSize(new Dimension(390, 124));
        jTextArea1.setEditable(false);
        jTextArea1.setMargin(new Insets(5, 5, 5, 5));
        jTextArea1.setText(wizardStep.getDescription());
        jTextArea1.setLineWrap(true);
        jTextArea1.setWrapStyleWord(true);

        yesBtn.setText("Yes");
        yesBtn.addActionListener(this);
        yesNoGroup.add(yesBtn);

        noBtn.setText("No");
        noBtn.addActionListener(this);
        yesNoGroup.add(noBtn);
        noBtn.setSelected(true);

        this.add(jTextArea1,  new GridBagConstraints(0, 0, 3, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(4, 5, 0, 5), 0, 12));
        this.add(noBtn,   new GridBagConstraints(2, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(10, 0, 13, 0), 20, 2));
        this.add(jLabel1,  new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(10, 24, 13, 0), 14, 9));
        this.add(yesBtn,   new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(10, 5, 13, 0), 14, 1));
        yesNoGroup.add(noBtn);
        yesNoGroup.add(yesBtn);
    }


    public void actionPerformed(ActionEvent e)
    {
        if(e.getSource() == yesBtn)
            wizard.enableNext();
        else if(e.getSource() == noBtn)
            wizard.disableNext();
    }

    public boolean getAnswer()
    {
        if(yesBtn.isSelected())
            return true;
        else
            return false;
    }
}
