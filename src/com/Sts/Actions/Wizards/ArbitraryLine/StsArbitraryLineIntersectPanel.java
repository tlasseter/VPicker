package com.Sts.Actions.Wizards.ArbitraryLine;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;

import javax.swing.*;
import java.awt.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsArbitraryLineIntersectPanel extends JPanel
{
    private StsArbitraryLineWizard wizard;
    private StsArbitraryLineIntersect wizardStep;

    StsGroupBox addBox = new StsGroupBox("Define Intersections");
    ButtonGroup orientGrp = new ButtonGroup();
    public StsToggleButton deleteIntersectionBtn;
    public StsToggleButton addIntersectionBtn;
    GridBagLayout gridBagLayout2 = new GridBagLayout();
    GridBagLayout gridBagLayout1 = new GridBagLayout();

    StsGroupBox intersectionBox = new StsGroupBox("Intersection Parameters");
    JLabel listLbl = new JLabel("Intersections:");
    JComboBox intersectionList = new JComboBox();
    JLabel rowLbl = new JLabel("InLine:");
    JTextField rowTxt = new JTextField();
    JLabel colLbl = new JLabel("CrossLine:");
    JTextField colTxt = new JTextField();

    static final String addButtonTip = "Create an intersecing line";
    static final String deleteButtonTip = "Select and delete an intersecting line";

    public StsArbitraryLineIntersectPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsArbitraryLineWizard)wizard;
        this.wizardStep = (StsArbitraryLineIntersect)wizardStep;

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
         this.setLayout(gridBagLayout2);
        addBox.setLayout(gridBagLayout1);

        addIntersectionBtn = new StsToggleButton("Add", addButtonTip, this, "addIntersection");
        deleteIntersectionBtn = new StsToggleButton("Delete", deleteButtonTip, this, "deleteIntersection");

        add(intersectionBox, new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0,
            GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 5, 0, 5), 0, 0));
        add(addBox,    new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(2, 5, 5, 5), 0, 0));

        addIntersectionBtn.setText("Add");
        addBox.add(addIntersectionBtn,   new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 0, 0, 0), 50, 10));
        deleteIntersectionBtn.setText("Delete");
        addBox.add(deleteIntersectionBtn,   new GridBagConstraints(1, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 0, 0, 0), 50, 10));


        intersectionBox.setLayout(new GridBagLayout());
        rowTxt.setEnabled(false);
        rowTxt.setBorder(BorderFactory.createEtchedBorder());
        rowTxt.setBackground(Color.lightGray);
        colTxt.setEnabled(false);
        colTxt.setBorder(BorderFactory.createEtchedBorder());
        colTxt.setBackground(Color.lightGray);

        intersectionBox.add(listLbl, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0,
           GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(5, 0, 0, 0), 0, 0));
        intersectionBox.add(intersectionList, new GridBagConstraints(1, 0, 1, 1, 1.0, 1.0,
           GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(5, 0, 0, 0), 100, 0));
        intersectionBox.add(rowLbl, new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0,
            GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(2, 0, 0, 0), 0, 0));
        intersectionBox.add(rowTxt, new GridBagConstraints(1, 1, 1, 1, 1.0, 1.0,
            GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 0, 0, 0), 100, 0));
        intersectionBox.add(colLbl, new GridBagConstraints(0, 2, 1, 1, 1.0, 1.0,
            GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(2, 0, 0, 0), 0, 0));
        intersectionBox.add(colTxt, new GridBagConstraints(1, 2, 1, 1, 1.0, 1.0,
            GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 0, 0, 0), 100, 0));

        this.setMinimumSize(new Dimension(300, 300));

    }

    public void initialize()
    {
        int nItems = initializeComboBox();
        addIntersectionBtn.setEnabled(true);
        deleteIntersectionBtn.setEnabled(nItems > 0);
    }

    private int initializeComboBox()
    {
        StsObjectRefList intersections = wizard.getIntersections();
        if(intersections == null) return 0;
        int nIntersections = intersections.getSize();
        if(nIntersections == 0)
        {
            ;
        }
        else
        {
            ;
        }
        return nIntersections;
    }

    public void deleteIntersection()
    {
        return;
    }

    public void addIntersection()
    {
        return;
    }

}
