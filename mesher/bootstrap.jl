(pwd() != @__DIR__) && cd(@__DIR__) # allow starting app from bin/ dir

using Mesher
const UserApp = Mesher
Mesher.main()
