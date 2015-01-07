package com.Sts.Actions.Wizards.CrossPlot;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;

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

public class StsSelectCrossplotPanel extends JPanel implements ActionListener
{
    private StsCrossplotWizard wizard;
    private StsSelectCrossplot wizardStep;

    private StsModel model = null;
    private StsCrossplot selectedCrossplot = null;

    JList crossplotList = new JList();
    DefaultListModel crossplotListModel = new DefaultListModel();

    JButton newCrossplotButton = new JButton();
    Border border1;

    StsCrossplot[] crossplots;
    GridBagLayout gridBagLayout1 = new GridBagLayout();

    public StsSelectCrossplotPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsCrossplotWizard)wizard;
        this.wizardStep = (StsSelectCrossplot)wizardStep;
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

    private void initialize()
    {
        model = wizard.getModel();
        crossplots = (StsCrossplot[])model.getCastObjectList(StsCrossplot.class);
        int nCrossplots = crossplots.length;
        for(int n = 0; n < nCrossplots; n++)
            crossplotListModel.addElement(crossplots[n].getName());
        crossplotList.setModel(crossplotListModel);
    }

    public StsCrossplot getSelectedCrossplot()
    {
        if(crossplotList.isSelectionEmpty()) return null;
        int selectedIndex = crossplotList.getSelectedIndex();
        return crossplots[selectedIndex];
    }

    void jbInit() throws Exception
    {
        border1 = BorderFactory.createCompoundBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED,Color.white,Color.white,new Color(178, 178, 178),new Color(124, 124, 124)),BorderFactory.createEmptyBorder(10,10,10,10));
        this.setLayout(gridBagLayout1);
        newCrossplotButton.setText("New Cross Plot");
        crossplotList.setBorder(BorderFactory.createEtchedBorder());
        crossplotList.setMaximumSize(new Dimension(200, 200));
        crossplotList.setMinimumSize(new Dimension(50, 50));
        crossplotList.setPreferredSize(new Dimension(200, 200));
        this.add(crossplotList,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 10, 0, 10), 0, 0));
        this.add(newCrossplotButton,   new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 0, 5, 0), 100, 0));
        newCrossplotButton.addActionListener(this);
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();

        // select new directory
        if(source == newCrossplotButton)
        {
            wizard.createNewCrossplot();
        }
    }
}
