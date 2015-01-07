package com.Sts.Actions.Wizards.Collaboration;

import com.Sts.Actions.*;
import com.Sts.DB.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.event.*;

public class StsSetCollaborationMode extends StsAction
{
    JDialog dialog;
    StsRadioButtonSelectPanel panel;
    StsDBFile db = null;

    static public final byte COLLABORATION_NONE = 0; // standard mode
    static public final byte COLLABORATION_MASTER = 1; // after loading database "db.*", write transactions to "db.*.n" where n is transaction number
    static public final byte COLLABORATION_SLAVE = 2; // after loading database "db.*", read transactions "db.*.n" where n is transaction number
    static public final String[] modeStrings = new String[] { "none", "master", "slave" };
    static public byte currentMode = COLLABORATION_NONE;
    static public CollaborationWindow collWindow;
    static public double bytesPerSec = 1250000; // 10Mbit LAN

    static public void setCurrentMode(byte mode) {

        currentMode = mode;

        // update the UI
        if(mode == COLLABORATION_MASTER) {
            collWindow.enableButton(false);
        } else {
            collWindow.enableButton(true);
        }
    }

    static public void addRequestListener(ActionListener al) {
        collWindow.addButtonListener(al);
    }

    public StsSetCollaborationMode(StsActionManager actionManager)
    {
        super(actionManager, true);
    }

    public boolean start()
    {
        try
        {
            db = model.getDatabase();
            /*
            if(db == null)
                new StsMessage(model.win3d, StsMessage.ERROR, "Must load a database before setting collaboration.");
            else */
            if (currentMode != COLLABORATION_NONE)
                new StsMessage(
                    model.win3d,
                    StsMessage.ERROR,
                    "Sorry, collaboration already set to: "
                        + modeStrings[currentMode]
                        + "\nCan't change: reload database.");
            else {

                StsCollaborationWizard collabWizard = new StsCollaborationWizard((StsActionManager) actionManager);
                collabWizard.start();

                dialog = new JDialog(model.win3d, "Set Collaboration mode", true);
                panel =
                    new StsRadioButtonSelectPanel(
                        "Set collaboration mode:",
                        new String[] { "None", "Master", "Slave" });
                panel.setButtonIndexSelected((int) currentMode);
                addButtonListeners();
                dialog.getContentPane().add(panel);
                dialog.pack();
                dialog.setVisible(true);
            }
            actionManager.endCurrentAction();
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsSetCollaborationMode.start() failed.",
                e, StsException.WARNING);
            return false;
        }
    }

    private void addButtonListeners()
    {
        ItemListener listener = new ItemListener()
        {
            public void itemStateChanged(ItemEvent event)
            {
		        byte buttonIndexSelected = (byte)panel.getButtonIndexSelected();
                if(currentMode != buttonIndexSelected) dialog.dispose();
                currentMode = buttonIndexSelected;
            }
        };
        panel.addButtonListeners(listener);
    }

    public boolean end()
    {
        try
        {
            /*
            StsSetCollaborationMode.collWindow = new CollaborationWindow();
            StsSetCollaborationMode.collWindow.setVisible(true);

            // start presence
            System.out.println("********* Starting Presence");
            BasicConfigurator.configure();
            Kernel k = null;

            KernelConfigContainer kcc = generateKernelConfigContainer();
            try {
                k = CommandReceiver.initializeKernel(kcc);
            } catch (AlreadyInitializedException aie) {
                System.err.println("Kernel already initialized");
                aie.printStackTrace();
            }

            System.out.println("******** Presence started");

            // if we are a master, create a pawn to handle all interactions
            if (currentMode == COLLABORATION_MASTER) {
                System.out.println("******** Creating DBPawn");
                DBPawn pawn = new DBPawn(db, model);

                pawn.regRoot(pawn);

                // add a hook to notify the pawn of changes in the model
                System.out.println("******** Inserting transaction hook");
                db.addTransactionListener(new PawnWriter(pawn));
            }

            // if we are a slave, connect appropriately
            if(currentMode == COLLABORATION_SLAVE) {
                DBPawn.setDB(db);
                DBPawn.setModel(model);
               DBPawn.setGLPanel(glPanel);
            }

            // show the SessionConnScheme
            SessionCSConfigContainer scscc = generateSCSCC("" + k.getUserID().getUniqueID());
            SessionConnectionScheme scc = null;

            try {
                scc = new SessionConnectionScheme(scscc, k);
                scc.setConnectionSchemeUI(new PlaceholderUI());
            } catch(Exception e) {
                System.err.println("Unable to create SessionConnectionScheme");
                e.printStackTrace();
            }

            try {
                if (currentMode == COLLABORATION_MASTER) {
                    scc.createSession("S2S", "");
                } else {
                    scc.joinSession(scc.getSessionByName("S2S"), "");
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
*/
            /*
            if(db == null) return true;
            if(currentMode == COLLABORATION_SLAVE) db.loadTransactions(model);
            return true;
            */
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsSetCollaborationMode.end() failed.",
                e, StsException.WARNING);
            currentMode = COLLABORATION_NONE;
            return false;
        }
    }
/*
    private SessionCSConfigContainer generateSCSCC(String uid) {

        SessionCSConfigContainer scscc = new SessionCSConfigContainer();
        System.out.println("UID = " + uid);
        scscc.setUniqueName(uid);
        scscc.setSessionType("S2S_SESSION");

        ClientConfigContainer ccc = new ClientConfigContainer();
        SettingsConfigContainer scc = new SettingsConfigContainer("default","MAUI",4353,"c:/collaboration/s2sdemo/infoserver/ar/infoserver/data/keystore","user","password");
        scc.setKeyStorePassphrase("passphrase");
        scc.setStrongSecurity(false);
        scc.setID(ClientConfigContainer.CURRENT_ID);

        ccc.add(scc);

        scscc.add(ccc);

        try {
            scscc.validate();
            return scscc;
        } catch(ValidationException ve) {
            ve.printStackTrace();
            return null;
        }
    }
*/
    /**
     * Method generateKernelConfigContainer.
     * @return KernelConfigContainer
     */
/*
    private KernelConfigContainer generateKernelConfigContainer() {

        KernelConfigContainer kcc = new KernelConfigContainer();

        // add a UserID ---------------
        kcc.add(new UserConfigContainer());

        // create the nengine config ------------ just use defaults.
        NetworkEngineConfigContainer necc = new NetworkEngineConfigContainer();
        kcc.add(necc);

        // add the BasicConnectivityManager configuration to the config
        BCMConfigContainer bcm = new BCMConfigContainer();

        bcm.setID(NetworkEngineConfigContainer.CM_ID);
        bcm.setGuaranteedPort(4600);
        bcm.setUDPPort(4000);
        bcm.setNegotiatorClassName("ar.connscheme.CSNegotiator");

        necc.setConnectivityManagerClassName(
            "ar.connmanager.basic.BasicConnectivityManager");
        necc.add(bcm);

        // create the console configs ------------
        ConsoleConfigContainer ccc = new ConsoleConfigContainer();

        // create the primary console
        ccc.setConsole(new SimpleConsole()); //new TextConsole());
        ccc.setID(KernelConfigContainer.PRIMARY_CONSOLE_ID);

        kcc.add(ccc);

        // ensure everything is configured properly
        try {
            kcc.validate();
        } catch (Exception ve) {
            System.err.println("Invalid configuration");
            ve.printStackTrace();
            return null;
        }

        return kcc;
    }
 */
}


//class PlaceholderUI implements ConnectionSchemeUI, SessionUI
//{
    /**
     * @see ar.connscheme.ConnectionSchemeUI#dispose()
     */
//    public void dispose() {
//    }

    /**
     * @see ar.connscheme.ConnectionSchemeUI#setVisible(boolean)
     */
//    public void setVisible(boolean visible) {
//    }

    /**
     * @see ar.sessionconnectionscheme.SessionUI#addParticipant(Participant)
     */
//    public void addParticipant(Participant name) {
//    }

    /**
     * @see ar.sessionconnectionscheme.SessionUI#addSession(Session)
     */
//    public void addSession(Session s) {
//    }

    /**
     * @see ar.sessionconnectionscheme.SessionUI#connectionLost()
     */
//    public void connectionLost() {
//    }

    /**
     * @see ar.sessionconnectionscheme.SessionUI#joinedSession(Session)
     */
//    public void joinedSession(Session s) {
//    }

    /**
     * @see ar.sessionconnectionscheme.SessionUI#joinFailed(Session)
     */
//    public void joinFailed(Session s) {
//    }

    /**
     * @see ar.sessionconnectionscheme.SessionUI#leftSession()
     */
//    public void leftSession() {
//    }

    /**
     * @see ar.sessionconnectionscheme.SessionUI#removeParticipant(Participant)
     */
//    public void removeParticipant(Participant name) {
//    }

    /**
     * @see ar.sessionconnectionscheme.SessionUI#removeSession(Session)
     */
//    public void removeSession(Session s) {
//    }

    /**
     * @see ar.sessionconnectionscheme.SessionUI#updateParticipant(Participant)
     */
//    public void updateParticipant(Participant name) {
//    }
//}

class CollaborationWindow extends JFrame {

    private JButton reqButton;

    public CollaborationWindow() {
        super("Collaboration Window");

        initButton();

        this.setSize(200,100);
    }

    private void initButton() {
        JPanel panel = new JPanel();

        reqButton = new JButton("Request Master");

        panel.add(reqButton);

        if(StsSetCollaborationMode.currentMode == StsSetCollaborationMode.COLLABORATION_MASTER) {
            reqButton.setEnabled(false);
        }

        this.setContentPane(panel);
    }

    void addButtonListener(ActionListener al) {
        this.reqButton.addActionListener(al);
    }

    void enableButton(boolean state) {
        this.reqButton.setEnabled(state);
    }
}

//class SimpleConsole extends Console {

    /**
     * @see presence.base.Console#printError(String, String)
     */
//    public void printError(String srcObj, String error) {
//        System.err.println(error);
//    }

    /**
     * @see presence.base.Console#println(String)
     */
//    public void println(String line) {
//        System.out.println(line);
//    }
//
//}
