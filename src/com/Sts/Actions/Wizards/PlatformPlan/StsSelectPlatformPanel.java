package com.Sts.Actions.Wizards.PlatformPlan;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;

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

public class StsSelectPlatformPanel extends JPanel implements ActionListener
{
    private StsPlatformPlanWizard wizard;
    private StsSelectPlatform wizardStep;

    private StsModel model = null;
    private StsPlatform selectedPlatform = null;

    JList platformList = new JList();
    DefaultListModel platformListModel = new DefaultListModel();

    JButton newPlatformButton = new JButton();

    StsPlatform[] platforms;
    GridBagLayout gridBagLayout1 = new GridBagLayout();

    public StsSelectPlatformPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsPlatformPlanWizard)wizard;
        this.wizardStep = (StsSelectPlatform)wizardStep;


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
        model = wizard.getModel();
        platforms = (StsPlatform[])model.getCastObjectList(StsPlatform.class);
        int nPlatforms = platforms.length;
        platformListModel.removeAllElements();
        for(int n = 0; n < nPlatforms; n++)
            platformListModel.addElement(platforms[n].getName());
        platformList.setModel(platformListModel);
    }

    public StsPlatform getSelectedPlatform()
    {
        if(platformList.isSelectionEmpty())
            return null;
        int selectedIndex = platformList.getSelectedIndex();
        return platforms[selectedIndex];
    }


    void jbInit() throws Exception
    {
        this.setLayout(gridBagLayout1);
        newPlatformButton.setText("New Platform");
        platformList.setBorder(BorderFactory.createEtchedBorder());
        platformList.setMaximumSize(new Dimension(200, 200));
        platformList.setMinimumSize(new Dimension(50, 50));
        platformList.setPreferredSize(new Dimension(200, 200));
        this.add(platformList,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.VERTICAL, new Insets(5, 10, 0, 10), 200, 0));
        this.add(newPlatformButton,   new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 0, 5, 0), 100, 0));
        newPlatformButton.addActionListener(this);
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();
        if(source == newPlatformButton)
        {
            wizard.createPlatform();
        }
    }

}
