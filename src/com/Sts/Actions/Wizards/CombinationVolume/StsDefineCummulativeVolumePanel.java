package com.Sts.Actions.Wizards.CombinationVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.Icons.*;
import com.Sts.UI.*;

import javax.swing.*;
import javax.swing.event.*;
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

public class StsDefineCummulativeVolumePanel extends JPanel implements ActionListener, ListSelectionListener
{
    private StsWizard wizard;
    private StsDefineCummulativeVolume wizardStep;

    private JPanel selectionPanel = new JPanel();
    private JPanel transferPanel = new JPanel();
    private JButton addBtn = new JButton();
    private JButton addAllBtn = new JButton();
    private JButton removeBtn = new JButton();
    private JButton removeAllBtn = new JButton();

    private StsBoxSetSubVolume selectedSet = null;
    private int selectedIdx = 0;
    private DefaultListModel availableVolsListModel = new DefaultListModel();
    private StsAbstractFileSet availableFileSet;
    private StsAbstractFile[] availableFiles;
    private JScrollPane availableVolsScrollPane = new JScrollPane();
    private JList availableVolsList = new JList();

    private DefaultListModel selectedVolsListModel = new DefaultListModel();
    private JScrollPane selectedVolsScrollPane = new JScrollPane();
    private JList selectedVolsList = new JList();

    JPanel jPanel1 = new JPanel();

    static final int UP = 0;
    static final int DOWN = 1;

    GridBagLayout gridBagLayout1 = new GridBagLayout();
    GridBagLayout gridBagLayout3 = new GridBagLayout();

    private StsModel model = null;
    ButtonGroup modeGroup = new ButtonGroup();
    JRadioButton sumBtn = new JRadioButton();
    JRadioButton averageBtn = new JRadioButton();
    JRadioButton firstBtn = new JRadioButton();
    JRadioButton lastBtn = new JRadioButton();
    GridBagLayout gridBagLayout2 = new GridBagLayout();
  JButton upBtn = new JButton();
  JButton downBtn = new JButton();
  GridBagLayout gridBagLayout4 = new GridBagLayout();
  GridBagLayout gridBagLayout5 = new GridBagLayout();

    public StsDefineCummulativeVolumePanel(StsWizard wizard, StsDefineCummulativeVolume wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = wizardStep;
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
        StsBoxSetSubVolume[] boxSets = (StsBoxSetSubVolume[])model.getCastObjectList(StsBoxSetSubVolume.class);

        int nSets = boxSets.length;
        availableVolsListModel.removeAllElements();
        for(int i=0; i<nSets; i++)
            availableVolsListModel.addElement(boxSets[i].getName());

        StsBoxSetSubVolumeClass boxSetClass = (StsBoxSetSubVolumeClass)model.getStsClass(StsBoxSetSubVolume.class);
        boxSetClass.setIsVisible(false);
    }


    void jbInit() throws Exception
    {
        this.setLayout(gridBagLayout5);
        selectionPanel.setBorder(BorderFactory.createEtchedBorder());
        selectionPanel.setLayout(gridBagLayout4);
//        this.setMinimumSize(new Dimension(0, 0));
//        this.setPreferredSize(new Dimension(0, 0));
        removeAllBtn.setFont(new java.awt.Font("Dialog", 0, 11));
        removeAllBtn.setMaximumSize(new Dimension(100, 20));
        removeAllBtn.setMinimumSize(new Dimension(100, 20));
        removeAllBtn.setPreferredSize(new Dimension(100, 20));
        removeAllBtn.setMargin(new Insets(0, 0, 0, 0));
        removeAllBtn.setText("<< Remove All");
        removeAllBtn.addActionListener(this);
        removeBtn.setFont(new java.awt.Font("Dialog", 0, 11));
        removeBtn.setMaximumSize(new Dimension(100, 20));
        removeBtn.setMinimumSize(new Dimension(100, 20));
        removeBtn.setPreferredSize(new Dimension(100, 20));
        removeBtn.setMargin(new Insets(0, 0, 0, 0));
        removeBtn.setText("< Remove");
        removeBtn.addActionListener(this);
        addAllBtn.setFont(new java.awt.Font("Dialog", 0, 11));
        addAllBtn.setMaximumSize(new Dimension(100, 20));
        addAllBtn.setMinimumSize(new Dimension(100, 20));
        addAllBtn.setPreferredSize(new Dimension(100, 20));
        addAllBtn.setMargin(new Insets(0, 0, 0, 0));
        addAllBtn.setText("Add All >>");
        addAllBtn.addActionListener(this);
        transferPanel.setLayout(gridBagLayout1);
//        transferPanel.setMaximumSize(new Dimension(100, 200));
//        transferPanel.setMinimumSize(new Dimension(25, 140));
//        transferPanel.setPreferredSize(new Dimension(50, 150));
//        selectedVolsList.setMinimumSize(new Dimension(75, 200));
        selectedVolsList.setModel(selectedVolsListModel);
        availableVolsList.setModel(availableVolsListModel);
        selectedVolsList.addListSelectionListener(this);
        availableVolsList.addListSelectionListener(this);
        addBtn.setText("Add >");
        addBtn.setFont(new java.awt.Font("Dialog", 0, 11));
        addBtn.setMaximumSize(new Dimension(100, 20));
        addBtn.setMinimumSize(new Dimension(100, 20));
        addBtn.setPreferredSize(new Dimension(100, 20));
        addBtn.setMargin(new Insets(0, 0, 0, 0));
        addBtn.addActionListener(this);
        availableVolsScrollPane.setPreferredSize(new Dimension(150, 110));
        selectedVolsScrollPane.setPreferredSize(new Dimension(150, 110));
        jPanel1.setBorder(BorderFactory.createEtchedBorder());
        jPanel1.setLayout(gridBagLayout2);
        sumBtn.setText("Sum");
        averageBtn.setSelected(true);
        averageBtn.setText("Average");
        firstBtn.setText("First");
        lastBtn.setText("Last");

        ImageIcon icon = StsIcon.createIcon("upArrow.gif");
        upBtn.setIcon(icon);
        icon = StsIcon.createIcon("downArrow.gif");
        downBtn.setIcon(icon);

        upBtn.addActionListener(this);
        downBtn.addActionListener(this);

        this.add(selectionPanel,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 4, 0, 5), 0, 2));
        transferPanel.add(removeBtn,        new GridBagConstraints(0, 2, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(20, 4, 0, 1), 0, 0));
        transferPanel.add(addAllBtn,      new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(20, 4, 0, 1), 0, 0));
        transferPanel.add(addBtn,      new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(4, 4, 0, 1), 0, 0));
        transferPanel.add(removeAllBtn,         new GridBagConstraints(0, 3, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(20, 4, 0, 1), 0, 0));
        selectionPanel.add(downBtn,  new GridBagConstraints(3, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 2, 5), 0, 0));
        selectionPanel.add(upBtn,  new GridBagConstraints(2, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 2, 0), 0, 0));

        this.add(jPanel1,  new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 4, 4, 5), 0, 9));

        selectionPanel.add(availableVolsScrollPane,  new GridBagConstraints(0, 0, 1, 2, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 5, 3, 0), 100, 127));
        selectionPanel.add(selectedVolsScrollPane,  new GridBagConstraints(2, 0, 2, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 0, 0, 4), 100, 98));

        selectedVolsScrollPane.getViewport().add(selectedVolsList, null);
        availableVolsScrollPane.getViewport().add(availableVolsList, null);

        selectionPanel.add(transferPanel,  new GridBagConstraints(1, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 0, 0, 0), 74, 60));

        jPanel1.add(sumBtn,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(2, 2, 4, 0), 0, 0));
        jPanel1.add(lastBtn,  new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(2, 48, 4, 13), 10, -3));
        jPanel1.add(averageBtn,  new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(2, 49, 4, 0), 2, -1));
        jPanel1.add(firstBtn,  new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(2, 49, 4, 0), 6, -5));

        modeGroup.add(averageBtn);
        modeGroup.add(sumBtn);
        modeGroup.add(firstBtn);
        modeGroup.add(lastBtn);
    }

    public void actionPerformed(ActionEvent e)
    {
        int i;
        int[] selectedIndices = null;

        Object source = e.getSource();

        selectedVolsList.removeListSelectionListener(this);

        // Add selected volumes
        if(source == addBtn)
        {
            selectedIndices = availableVolsList.getSelectedIndices();
            for(i=0; i< selectedIndices.length; i++)
            {
                if(selectedVolsListModel.indexOf(availableVolsList.getModel().getElementAt(selectedIndices[i])) >= 0)
                    continue;
                selectedVolsListModel.addElement(availableVolsList.getModel().getElementAt(selectedIndices[i]));
            }
        }
        // Remove selected volume(s)
        else if (source == removeBtn)
        {
            selectedIndices = selectedVolsList.getSelectedIndices();
            for(i= selectedIndices.length - 1; i>=0; i--)
                selectedVolsListModel.removeElementAt(selectedIndices[i]);
        }
        // Remove all volumes
        else if (source == removeAllBtn)
            selectedVolsListModel.removeAllElements();

        // Add all volumes
        else if (source == addAllBtn)
        {
            selectedVolsListModel.removeAllElements();
            for(i = 0; i < availableVolsList.getModel().getSize(); i++)
                selectedVolsListModel.addElement(availableVolsList.getModel().getElementAt(i));
        }
        // Up Button
        else if (source == upBtn)
        {
            if(selectedSet == null)
                new StsMessage(model.win3d, StsMessage.WARNING,"Must select item to move in left list.");
            else
            {
                moveSelectedItem(UP);
            }
        }
        // Down Button
        else if (source == downBtn)
        {
            if(selectedSet == null)
                new StsMessage(model.win3d, StsMessage.WARNING, "Must select item to move in left list.");
            else
            {
                moveSelectedItem(DOWN);
            }
        }
        selectedVolsList.addListSelectionListener(this);
    }

    public void valueChanged(ListSelectionEvent e)
    {
        Object source = e.getSource();
        int[] selectedIndices = null;
        String volName = null;

        if(source == selectedVolsList)
        {
            selectedIndices = selectedVolsList.getSelectedIndices();
            if(selectedIndices != null)
                volName = (String) selectedVolsList.getModel().getElementAt(selectedIndices[selectedIndices.length - 1]);
            else
                volName = null;
        }
        else if(source == availableVolsList)
        {
            selectedIndices = availableVolsList.getSelectedIndices();
            volName = (String)availableVolsList.getModel().getElementAt(selectedIndices[selectedIndices.length-1]);
        }
        selectedIdx = selectedIndices[selectedIndices.length-1];
        selectedSet = (StsBoxSetSubVolume) model.getObjectWithName(StsBoxSetSubVolume.class, volName);

        return;
    }

    public int getMode()
    {
        if(sumBtn.isSelected())
            return wizardStep.SUM;
        else if(averageBtn.isSelected())
            return wizardStep.AVERAGE;
        else if(firstBtn.isSelected())
            return wizardStep.FIRST;
        else if(lastBtn.isSelected())
            return wizardStep.LAST;
        else
            return -1;
    }

    private void moveSelectedItem(int dir)
    {
        StsBoxSetSubVolume[] list = null;
        int i;

        list = getSetList();
        if(list == null)
            return;
        switch(dir)
        {
            case UP:
                if(selectedIdx != 0)
                {
                    selectedVolsListModel.removeAllElements();
                    for (i = 0; i < selectedIdx-1; i++)
                        selectedVolsListModel.addElement(list[i].getName());
                    selectedVolsListModel.addElement(list[selectedIdx].getName());
                    selectedVolsListModel.addElement(list[selectedIdx-1].getName());
                    for (i = selectedIdx+1; i < list.length; i++)
                        selectedVolsListModel.addElement(list[i].getName());
                    setSelected(selectedIdx - 1);
                }
                break;
            case DOWN:
                if(selectedIdx != (list.length-1))
                {
                    selectedVolsListModel.removeAllElements();
                    for (i = 0; i < selectedIdx; i++)
                        selectedVolsListModel.addElement(list[i].getName());
                    selectedVolsListModel.addElement(list[selectedIdx+1].getName());
                    selectedVolsListModel.addElement(list[selectedIdx].getName());
                    for (i = selectedIdx+2; i < list.length; i++)
                        selectedVolsListModel.addElement(list[i].getName());
                    setSelected(selectedIdx + 1);
                }

                break;
        }
    }

    private void setSelected(int idx)
    {
        selectedIdx = idx;
        selectedSet = (StsBoxSetSubVolume) model.getObjectWithName(StsBoxSetSubVolume.class, (String) selectedVolsList.getModel().getElementAt(selectedIdx));
        selectedVolsList.setSelectedIndex(selectedIdx);
    }

    public StsBoxSetSubVolume[] getSetList()
    {
        String setName = null;
        StsBoxSetSubVolume[] boxSets = new StsBoxSetSubVolume[selectedVolsList.getModel().getSize()];
        for(int i=0; i<selectedVolsList.getModel().getSize(); i++)
        {
            setName = (String)selectedVolsList.getModel().getElementAt(i);
            boxSets[i] = (StsBoxSetSubVolume) model.getObjectWithName(StsBoxSetSubVolume.class, setName);
        }
        return boxSets;
    }
}
