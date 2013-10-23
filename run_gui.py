import OpenPNM

import flow_gui

# flow_gui.splash('splash_loading.png')
flow_gui.add_source(OpenPNM.GUI.network)
flow_gui.add_source(OpenPNM.Physics.CapillaryPressure, name="physics")
flow_gui.add_source(OpenPNM.GUI.algorithm)
flow_gui.add_source(OpenPNM.GUI.output)
flow_gui.set_previews(OpenPNM.GUI.plots)
flow_gui.run()



    # # Create and display the splash screen
    # splash_pix = QPixmap('splash_loading.png')
    # splash = QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)
    # splash.setMask(splash_pix.mask())
    # splash.show()
    # app.processEvents()

    # # Simulate something that takes time
    # time.sleep(2)

    # form = Form()
    # form.show()
    # splash.finish(form)
    # app.exec_()