package com.Sts.Actions.Wizards.Color;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
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
 * @author not attributable
 * @version 1.0
 */

public class StsSpectrumSelect extends JDialog implements ActionListener, ListSelectionListener
{
    private StsModel model = null;
    private StsSpectrum spectrum = null;
    private StsColorscalePanel csPanel = null;
    private StsColorscalePanel origPanel = null;
    private DefaultListModel availableSpectrumListModel = new DefaultListModel();
    private boolean success = false;

    JPanel jPanel1 = new JPanel();
    JPanel jPanel2 = new JPanel();
    JPanel selectionPanel = new JPanel();
    JList availableSpectrumList = new JList();
    JButton okBtn = new JButton();
    JScrollPane availablePalettesScrollPane = new JScrollPane();
    GridBagLayout gridBagLayout3 = new GridBagLayout();
    JPanel paletteViewPanel = new JPanel();
    JButton cancelBtn = new JButton();
    GridBagLayout gridBagLayout1 = new GridBagLayout();
  GridBagLayout gridBagLayout2 = new GridBagLayout();
  GridBagLayout gridBagLayout4 = new GridBagLayout();

    public StsSpectrumSelect(StsColorscalePanel panel)
    {
        super((Frame)null, "Select Color Palette", true);

        model = StsObject.getCurrentModel();
        origPanel = panel;

        try
        {
            initialize();
            jbInit();
            pack();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }


    private void initialize()
    {
        setAvailableSpectrums();
        availableSpectrumList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        String name = availableSpectrumListModel.getElementAt(0).toString();
        spectrum = model.getSpectrum(name);
        if(csPanel != null)
        {
            csPanel.setSpectrum(spectrum);
        }
        else
        {
            csPanel = new StsColorscalePanel(spectrum);
        }
    }

    private void setAvailableSpectrums()
    {
        availableSpectrumListModel.removeAllElements();
        StsClass spectrumClass = model.getCreateStsClass(StsSpectrum.class);
        int nSpectrums = spectrumClass.getSize();
        for(int n=0; n<nSpectrums; n++)
        {
            StsSpectrum spectrum = (StsSpectrum)spectrumClass.getElement(n);
            if(spectrum.getNColors() == 255)
                availableSpectrumListModel.addElement(spectrum.getName());
        }
    }

    private void jbInit() throws Exception
    {
        jPanel1.setLayout(gridBagLayout4);
        jPanel2.setLayout(gridBagLayout1);
        selectionPanel.setBorder(BorderFactory.createEtchedBorder());
        selectionPanel.setLayout(gridBagLayout2);
        availableSpectrumList.addListSelectionListener(this);
        okBtn.setFont(new java.awt.Font("Dialog", 0, 11));
        okBtn.setText("Ok");
        okBtn.addActionListener(this);
        availablePalettesScrollPane.setPreferredSize(new Dimension(150, 110));
        paletteViewPanel.setLayout(gridBagLayout3);
        cancelBtn.setText("Cancel");
        cancelBtn.addActionListener(this);
        availableSpectrumList.setModel(availableSpectrumListModel);
        paletteViewPanel.setMaximumSize(new Dimension(30, 449));
        paletteViewPanel.setMinimumSize(new Dimension(30, 449));
        paletteViewPanel.setPreferredSize(new Dimension(30, 300));
        paletteViewPanel.setToolTipText("");
        jPanel2.setBorder(BorderFactory.createEtchedBorder());
        this.getContentPane().add(jPanel1, BorderLayout.CENTER);
        selectionPanel.add(availablePalettesScrollPane,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(4, 3, 2, 0), -11, 194));
    selectionPanel.add(paletteViewPanel,  new GridBagConstraints(1, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(4, 18, 2, 5), 0, -4));
    paletteViewPanel.add(csPanel);
        jPanel1.add(jPanel2,  new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 2, 4, 5), -4, -3));
        jPanel2.add(okBtn,   new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(4, 5, 6, 0), 33, 0));
        jPanel2.add(cancelBtn,   new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(4, 29, 6, 5), 13, 0));
        availablePalettesScrollPane.getViewport().add(availableSpectrumList, null);
        jPanel1.add(selectionPanel,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(3, 2, 0, 5), 0, 0));
    }

    public void valueChanged(ListSelectionEvent e)
    {
        Object source = e.getSource();
        String name = null;

        name = availableSpectrumListModel.getElementAt(availableSpectrumList.getSelectedIndex()).toString();
        spectrum = model.getSpectrum(name);

        csPanel.setSpectrum(spectrum);
    }

    public void actionPerformed(ActionEvent e)
    {
        int i;
        Object source = e.getSource();

        if (source == okBtn)
        {
            int selectedIndex = availableSpectrumList.getSelectedIndex();
			if(selectedIndex < 0) return;
            spectrum = model.getSpectrum(availableSpectrumList.getModel().getElementAt(selectedIndex).toString());
            if(origPanel != null) origPanel.setSpectrum(spectrum);
            this.setVisible(false);
            success = true;
        }
        else if (source == cancelBtn)
        {
            this.setVisible(false);
            success = false;
        }
    }

    public StsSpectrum getSelectedSpectrum()
    {
        return spectrum;
    }

    public boolean getSuccess()
    {
        return success;
    }
}
