
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Super;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

import java.awt.*;

public class StsSuperGeneral extends StsWizardStep
{
    StsSuperGeneralPanel panel;
    StsHeaderPanel header;
    String question = "Not Set";
    String description = "Not Set";
    int type = -1;

    public StsSuperGeneral(StsWizard wizard, int type, String title, String question, String description)
    {
        super(wizard);
        this.type = type;
        this.question = question;
        this.description = description;

        panel = new StsSuperGeneralPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(400, 400));
        header.setTitle(title);
        header.setSubtitle("Questions & Answers");
    }

    public boolean start()
    {
        return true;
    }

    public boolean end()
    {
        return true;
    }

    public boolean getAnswer()
    {
        return panel.getAnswer();
    }

    public int getType() { return type; }
    public String getQuestion() { return question; }
    public String getDescription() { return description; }
}

