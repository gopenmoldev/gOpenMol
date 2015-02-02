namespace eval BallPlugin {

;
# ::gom::EnablePlugin will load the plugin (by
# calling "::gom::LoadPlugin $plugin") and call
# a procedure "::[string totitle $plugin]Plugin::Enable" if it exists.
proc Enable {} {
    puts "A Tcl procedure '[namespace current]::Enable' called!"
}

# This will create a menu item "Enable ball plugin".  After a user
# has selected the menu item, a message "Ball plugin is enabled" is
# displayed.
gom::RegisterPlugin "ball" "ball plugin"

}
