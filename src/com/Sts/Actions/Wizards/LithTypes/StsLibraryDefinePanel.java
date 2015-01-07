package com.Sts.Actions.Wizards.LithTypes;

import com.Sts.Actions.Wizards.*;

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

public class StsLibraryDefinePanel extends JPanel implements ActionListener
{
    private StsLithTypesWizard wizard;
    private StsLibraryDefine wizardStep;

    JPanel jPanel1 = new JPanel();
    JLabel jLabel1 = new JLabel();
    JTextField libraryNameTxt = new JTextField();
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    GridBagLayout gridBagLayout2 = new GridBagLayout();

    public StsLibraryDefinePanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsLithTypesWizard)wizard;
        this.wizardStep = (StsLibraryDefine)wizardStep;
        try
        {
            jbInit();
            initialize();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initialize()
    {

    }

    void jbInit() throws Exception
    {
        this.setLayout(gridBagLayout2);
        jPanel1.setBorder(BorderFactory.createEtchedBorder());
        jPanel1.setLayout(gridBagLayout1);
        jLabel1.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel1.setText("Library Name:");
        libraryNameTxt.setText("libraryName");
        this.add(jPanel1,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 5, 8, 5), 3, 0));
        jPanel1.add(libraryNameTxt,  new GridBagConstraints(1, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(8, 0, 12, 5), 213, 4));
        jPanel1.add(jLabel1,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(8, 7, 12, 0), 28, 6));

    }

    public void actionPerformed(ActionEvent e)
    {

    }

    public String getLibraryName()
    {
        return libraryNameTxt.getText();
    }
}
