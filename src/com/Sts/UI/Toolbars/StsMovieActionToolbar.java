
/**
 * <p>Title: S2S Development</p>
 * <p>Description: Movie toolbar class. Defines the toolbar that is presented
 * once movies have been defined via the movie workflow step.</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author Stuat A. Jackson
 * @version 1.0
 */

package com.Sts.UI.Toolbars;

import com.Sts.Actions.Movie.*;
import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.border.*;
import java.awt.*;
import java.awt.event.*;

public class StsMovieActionToolbar extends StsToolbar implements StsSerializable
{
    /** Toolbar name used to reference the toolbar throughout the code */
    public static final String NAME = "Movie Action Toolbar";
    public static final boolean defaultFloatable = true;

    /** Unique identifier button name to start the movie - moviePlay */
    public static final String PLAY_ACTION = "moviePlay";
    public static final int PLAY = 0;
    /** Unique identifier button name to stop the movie and return to start - movieStop */
    public static final String STOP_ACTION = "movieStop";
    public static final int STOP = 1;
    /** Unique identifier button name to move to beginning of movie - movieStart */
    public static final String START_ACTION = "movieStart";
    public static final int START = 2;
    /** Unique identifier button name to move to end of movie - movieEnd */
    public static final String END_ACTION = "movieEnd";
    public static final int END = 3;
    /** Unique identifier button name to pause the move - moviePause */
    public static final String PAUSE_ACTION = "moviePause";
    public static final int PAUSE = 4;
    /** Unique identifier button name to play movie in reverse - movieReverese */
    public static final String REVERSE_ACTION = "movieReverse";
    public static final int REVERSE = 5;
    /** Unique identifier button name to save the movie frames as gif files - movieSave */
    public static final String SAVE_ACTION = "movieSave";
    public static final int SAVE = 6;

    /** The action object used to handle toolbar events */
    transient public StsMovieAction movieAction = null;
    transient public JPanel comboBoxPanel = null;
    transient public StsModel model = null;

    /** Developer Controlled Buttons */
    transient private StsToggleButton revBtn = null;
    transient private StsToggleButton playBtn = null;
    transient private StsToggleButton pauseBtn = null;
    transient private StsToggleButton stopBtn = null;
    transient ButtonGroup btnGroup = new ButtonGroup();

    /**
     * Movie toolbar constructor
     */
    public StsMovieActionToolbar()
     {
         super(NAME);
     }

    public StsMovieActionToolbar(StsWin3dBase win3d, StsMovieAction movieAction)
    {
        super(NAME);
        this.movieAction = movieAction;
        initialize(win3d);
    }

    public boolean initialize(StsWin3dBase win3d)
    {
        this.model = win3d.getModel();
        if(movieAction == null)
        {
            StsClass stsClass = model.getCreateStsClass(StsMovie.class);
            StsMovie movie = (StsMovie)stsClass.getCurrentObject();
            movieAction = new StsMovieAction(model, movie);
        }
        Border border = BorderFactory.createEtchedBorder();
        setBorder(border);

        JComboBox comboBox = createComboBox();
        comboBoxPanel = new JPanel(new GridBagLayout());
        comboBoxPanel.add(comboBox,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));

        add(comboBoxPanel);
        addSeparator();
        add(new StsButton(START_ACTION, "Return to Start of Movie", this, "startAction"));
        revBtn = new StsToggleButton(REVERSE_ACTION, "Play Movie in Reverse", this, "reverseAction");
        add(revBtn);
        stopBtn = new StsToggleButton(STOP_ACTION, "Stop and Return to Start of Movie", this, "stopAction");
        add(stopBtn);
        pauseBtn = new StsToggleButton(PAUSE_ACTION, "Pause on Current Frame", this, "pauseAction");
        add(pauseBtn);
        playBtn = new StsToggleButton(PLAY_ACTION, "Play from Current Frame", this , "playAction");
        add(playBtn);
        add(new StsButton(END_ACTION, "Go to the End of the Movie", this, "endAction"));
//        add(new StsButton(SAVE_ACTION, "Save a Set of Frames for the Movie", movieAction, "saveAction"));

        btnGroup.add(playBtn);
        btnGroup.add(revBtn);
        btnGroup.add(pauseBtn);
        btnGroup.add(stopBtn);

        addSeparator();
        addCloseIcon(model.win3d);

        setMinimumSize();
        return true;
    }

    public void startAction()
    {
        if(!initializeMovieAction()) return;
        movieAction.startAction();
    }

    public void reverseAction()
    {
        if(!initializeMovieAction()) return;
        movieAction.reverseAction();
    }

    public void stopAction()
    {
        if(!initializeMovieAction()) return;
        movieAction.stopAction();
    }

    public void pauseAction()
    {
        if(!initializeMovieAction()) return;
        movieAction.pauseAction();
    }

    public void playAction()
    {
        if(!initializeMovieAction()) return;
        movieAction.playAction();
    }

    public void endAction()
    {
        if(!initializeMovieAction()) return;
        movieAction.endAction();
    }

    private boolean initializeMovieAction()
    {
        if(movieAction != null) return true;
        StsClass stsClass = model.getCreateStsClass(StsMovie.class);
        StsMovie movie = (StsMovie)stsClass.getCurrentObject();
        if(movie == null) return false;
        movieAction = new StsMovieAction(model, movie);
        return true;
    }

    public StsMovieAction getMovieAction() { return movieAction; }
    private JComboBox createComboBox()
    {
        StsClass stsClass = model.getCreateStsClass(StsMovie.class);
        StsObject[] objects = stsClass.getObjectList();
        int nObjects = objects.length;
        DefaultComboBoxModel comboBoxModel = new DefaultComboBoxModel(objects);
        JComboBox comboBox = new JComboBox(comboBoxModel);
        comboBox.setName(StsMovie.class.getName());
        comboBox.setVisible(nObjects > 0);
        comboBox.setLightWeightPopupEnabled(false);

        StsObject currentObject = stsClass.getCurrentObject();
        if(nObjects > 0)
        {
            if(currentObject != null)
                comboBox.setSelectedItem(currentObject);
            else
                comboBox.setSelectedIndex(nObjects - 1);
        }

        stsClass.addActionListener(new ActionListener()
        {
            public void actionPerformed(ActionEvent e)
            {
                String command = e.getActionCommand();
                Object object = e.getSource();
                if(command.equals("add")) comboBoxAddObject(object);
                else if(command.equals("delete")) comboBoxDeleteObject(object);
                else if(command.equals("selected")) comboBoxSetItem(object);
            }
        });

        comboBox.addItemListener(new ItemListener()
        {
            public void itemStateChanged(ItemEvent e)
            {
                int state = e.getStateChange();
                Object object = e.getItem();
                if(state == ItemEvent.SELECTED) comboBoxSelectObject(object);
            }
        });
        return comboBox;
    }

    /*
     *  Called when comboBox item has been selected. Propagate this change
     *  to other interested parties and repaint.
     */
    private void comboBoxSelectObject(Object object)
    {
        JComboBox comboBox = getComboBox(object);
        movieAction.setMovie((StsMovie)object);
        if(comboBox == null) return;
        objectPanelSelectObject(object);
        StsObject stsObject = (StsObject)object;
        StsClass stsClass = model.getCreateStsClass(stsObject);
        stsClass.setCurrentObject(stsObject);
    }

    /**
     *  Called when a comboBox item has been changed externally and we need to
     *  change the item displayed in the comboBox.  By changing the model, we do
     *  not fire an item stateChanged causing an endless loop.
     */
    public void comboBoxSetItem(Object object)
    {
        JComboBox comboBox = getComboBox(object);
        if(comboBox == null) return;
        DefaultComboBoxModel comboBoxModel = (DefaultComboBoxModel)comboBox.getModel();
        comboBoxModel.setSelectedItem(object);
    }

    // if added, only update the view in the main window
    private void comboBoxAddObject(Object object)
    {
        if(comboBoxPanel == null) return;
        JComboBox comboBox = (JComboBox)comboBoxPanel.getComponent(0);
        if(comboBox == null) return;

        DefaultComboBoxModel comboBoxModel = (DefaultComboBoxModel)comboBox.getModel();

        comboBoxModel.addElement(object);
        comboBoxModel.setSelectedItem(object);
        if(comboBox.isVisible()) return;
        comboBox.setVisible(true);
        comboBoxPanel.setVisible(true);
        if(!isVisible()) setVisible(true);
    }

    private JComboBox getComboBox(Object object)
    {
        if(comboBoxPanel == null) return null;
        return (JComboBox)comboBoxPanel.getComponent(0);
    }
    private void comboBoxDeleteObject(Object object)
    {
        if(comboBoxPanel == null) return;
        JComboBox comboBox = (JComboBox)comboBoxPanel.getComponent(0);
        if(comboBox == null) return;

        DefaultComboBoxModel comboBoxModel = (DefaultComboBoxModel)comboBox.getModel();
        comboBoxModel.removeElement(object);
        comboBox.setSelectedIndex(0);
        object = comboBox.getSelectedItem();
    }

    private void objectPanelSelectObject(Object object)
    {
        if(model.win3d.isMainWindow()) model.win3d.objectTreePanel.selected((StsObject)object);
    }

    public void setState(int state_)
    {
        final int state = state_;
        StsToolkit.runLaterOnEventThread
        (
            new Runnable()
            {
                public void run()
                {
                    doSetState(state);
                }
            }
        );
    }

    public void doSetState(int state)
    {
        switch(state)
        {
            case START:
            case END:
                stopBtn.setSelected(true);
                revBtn.setEnabled(true);
                playBtn.setEnabled(true);
                break;
            case STOP:
                stopBtn.setSelected(true);
                revBtn.setEnabled(true);
                playBtn.setEnabled(true);
                break;
            case PAUSE:
                revBtn.setEnabled(true);
                playBtn.setEnabled(true);
                break;
            case PLAY:
                playBtn.setSelected(true);
                revBtn.setEnabled(false);
                break;
            case REVERSE:
                revBtn.setSelected(true);
                playBtn.setEnabled(false);
                break;
        }
        repaint();
    }
}

