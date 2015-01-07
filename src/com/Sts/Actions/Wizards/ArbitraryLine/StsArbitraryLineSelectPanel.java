package com.Sts.Actions.Wizards.ArbitraryLine;

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

public class StsArbitraryLineSelectPanel extends JPanel implements ActionListener
{
    private StsArbitraryLineWizard wizard;
    private StsArbitraryLineSelect wizardStep;

    private StsModel model = null;
    private StsArbitraryLine selectedLine = null;

    JList lineList = new JList();
    DefaultListModel lineListModel = new DefaultListModel();

    JButton newLineButton = new JButton();
    Border border1;

    StsArbitraryLine[] lines;
  JLabel jLabel1 = new JLabel();
  GridBagLayout gridBagLayout1 = new GridBagLayout();

    public StsArbitraryLineSelectPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsArbitraryLineWizard)wizard;
        this.wizardStep = (StsArbitraryLineSelect)wizardStep;
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
        lines = (StsArbitraryLine[])model.getCastObjectList(StsArbitraryLine.class);
        int nLines = lines.length;
        for(int n = 0; n < nLines; n++)
            lineListModel.addElement(lines[n].getName());
        lineList.setModel(lineListModel);
    }

    public StsArbitraryLine getSelectedLine()
    {
        if(lineList.isSelectionEmpty()) return null;
        int selectedIndex = lineList.getSelectedIndex();
        return lines[selectedIndex];
    }

    void jbInit() throws Exception
    {
        border1 = BorderFactory.createCompoundBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED,Color.white,Color.white,new Color(178, 178, 178),new Color(124, 124, 124)),BorderFactory.createEmptyBorder(10,10,10,10));
        this.setLayout(gridBagLayout1);
        newLineButton.setText("New Line");
        lineList.setBorder(BorderFactory.createEtchedBorder());
        lineList.setMaximumSize(new Dimension(200, 200));
        lineList.setMinimumSize(new Dimension(50, 50));
        lineList.setPreferredSize(new Dimension(200, 200));
        jLabel1.setText("Available Lines");
        this.add(lineList,  new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 100, 0, 100), 0, 0));
        this.add(newLineButton,      new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(38, 140, 4, 140), 16, 0));
        this.add(jLabel1,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(8, 100, 0, 100), 126, 9));
        newLineButton.addActionListener(this);
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();

        // select new directory
        if(source == newLineButton)
        {
            wizard.createNewLine();
        }
    }
}
