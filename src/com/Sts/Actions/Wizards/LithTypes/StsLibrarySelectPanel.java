package com.Sts.Actions.Wizards.LithTypes;

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

public class StsLibrarySelectPanel extends JPanel implements ActionListener
{
    private StsLithTypesWizard wizard;
    private StsLibrarySelect wizardStep;

    private StsModel model = null;
    private StsTypeLibrary selectedLibrary = null;

    JList libraryList = new JList();
    DefaultListModel libraryListModel = new DefaultListModel();

    JButton newLibraryButton = new JButton();
    Border border1;
    GridBagLayout gridBagLayout1 = new GridBagLayout();

    StsTypeLibrary[] libraries;

    public StsLibrarySelectPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsLithTypesWizard)wizard;
        this.wizardStep = (StsLibrarySelect)wizardStep;
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
        checkCreateDefaultLibraries();
        libraries = (StsTypeLibrary[])model.getCastObjectList(StsTypeLibrary.class);
        int nLibraries = libraries.length;
        if(nLibraries > 0)
        {
            for (int n = 0; n < nLibraries; n++)
                libraryListModel.addElement(libraries[n].getName());
            libraryList.setModel(libraryListModel);
            wizard.setSelectedLibrary(libraries[0]);
        }
    }

    private void checkCreateDefaultLibraries()
    {
        StsTypeLibrary.getCreateDefaultLibrary();
        StsTypeLibrary.getCreateGenericLibrary();
    }

    public void reinitialize()
    {
        libraryListModel.removeAllElements();

        model = wizard.getModel();
        libraries = (StsTypeLibrary[])model.getCastObjectList(StsTypeLibrary.class);
        int nLibraries = libraries.length;
        for(int n = 0; n < nLibraries; n++)
            libraryListModel.addElement(libraries[n].getName());
        libraryList.setModel(libraryListModel);
    }

    public StsTypeLibrary getSelectedLibrary()
    {
        if(libraryList.isSelectionEmpty()) return null;
        int selectedIndex = libraryList.getSelectedIndex();
        return libraries[selectedIndex];
    }

    void jbInit() throws Exception
    {
        border1 = BorderFactory.createCompoundBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED,Color.white,Color.white,new Color(178, 178, 178),new Color(124, 124, 124)),BorderFactory.createEmptyBorder(10,10,10,10));
        this.setLayout(gridBagLayout1);
        newLibraryButton.setText("New Library");
        libraryList.setBorder(BorderFactory.createEtchedBorder());
        libraryList.setMaximumSize(new Dimension(200, 200));
        libraryList.setMinimumSize(new Dimension(50, 50));
        libraryList.setPreferredSize(new Dimension(200, 200));
        this.add(libraryList,   new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(4, 4, 4, 4), 0, 0));
        this.add(newLibraryButton,   new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(4, 4, 4, 4), 0, 0));
        newLibraryButton.addActionListener(this);
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();

        // select new directory
        if(source == newLibraryButton)
        {
            wizard.createNewLibrary();
        }
    }
}
