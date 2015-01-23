#! /usr/bin/python3

### Importing Python modules

# Core modules
from functools import reduce
import glob, itertools, multiprocessing, operator, os, platform, re, shutil, string, subprocess, sys

# External modules
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from rpy2.robjects import r
from rpy2.robjects.packages import importr
from scipy.integrate import quad
import numpy as np
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri
import psutil

### All preliminary subroutines

## Calculate average length routine
def AveLength(sequencefile):

	# Starting variables
	count = 0
	total = 0

	# Processing fasta file
	origfasta = open(sequencefile, "rU") if os.path.isfile(sequencefile) else sys.exit("There is no way to open %s\n" % sequencefile)
	for record in SeqIO.parse(origfasta, "fasta"):
		count = count + SeqIO.write(record, "/dev/null", "fasta")
		total = total + len(record)
	origfasta.close()

	# Writing results
	average = int(total / count)
	return average

## Joining files (Unix command: "Cat") function
def Cat(joinedfile, listfile):
	alloldfiles = listfile.split(",")
	with open(joinedfile, "w") as out:
		for newfile in alloldfiles:
			with open(newfile) as single:
				out.write(single.read())
				
## Counting number of Headers routine
def CountingHeaders(fastafile):

	# Starting variables
	nrheaders = 0
	
	# Processing fasta file
	origfasta = open(fastafile, "rU") if os.path.isfile(fastafile) else sys.exit("There is no way to open %s\n" % fastafile)
	for record in SeqIO.parse(origfasta, "fasta"):
		nrheaders = nrheaders + SeqIO.write(record, "/dev/null", "fasta")
	origfasta.close()

	# Writing results	
	return nrheaders

## Declustering Fasta files routine
def DeclusterFASTAfile(fastafile, namefile, derepfile):

	newffile = fastafile.replace('.fasta','.declustered.fasta').replace('.clust','').replace('.derep','').replace('.sorted','')

	# Starting new dictionaries
	abundanceofcluster = {}
	latestsequencesheaders = {}
	nameofcluster = {}
	nameSeq = {}
	newsequencesheaders = {}
	nrseqheaders = {}
	oldsequencesheaders = {}
	sample = {}

	# Saving the DEREP file in memory
	my_derepfile = open(derepfile) if os.path.isfile(derepfile) else sys.exit("%s can not be found\n" % derepfile)
	derepline = my_derepfile.read().splitlines()
	for line in derepline:
		mainderepline = re.split(r"\t", line)
		refName = mainderepline[0]
		Synonims = mainderepline[1]
		oldsequencesheaders[refName] = Synonims
	my_derepfile.close()

	# Saving the NAME file in memory
	my_namefile = open(namefile) if os.path.isfile(namefile) else sys.exit("%s can not be found\n" % namefile)
	nameline = my_namefile.read().splitlines()
	for line in nameline:
		mainnameline = re.split(r"\t", line)
		SRefName = re.sub(';size=\d+;','',mainnameline[0])
		anotherelements = re.sub(';size=\d+;','',mainnameline[1])
		newsequencesheaders[SRefName] = anotherelements
	my_namefile.close()

	# Synonim table
	for key in newsequencesheaders.keys():
		if re.search(r",", newsequencesheaders[key]):
			marray = newsequencesheaders[key].rsplit(",")
		else:
			marray = [newsequencesheaders[key]]
		latestsequencesheaders[key] = []
		for element in marray:
			element = element.replace(element, oldsequencesheaders[element])
			if re.search(r",", element):			
				templist = element.rsplit(",")
			else:
				templist = [element]
			templist = str(templist).replace("[", "").replace("]", "").replace("'","").replace(" ", "")
			if re.search(r",", templist):
				newarray = templist.rsplit(",") 
			else:
				newarray = [templist]
			latestsequencesheaders[key].extend(newarray)
			nrseqheaders[key] = len(newarray)

	# Analysing and printing the declustered FASTA file
	origfasta = open(fastafile, "rU") if os.path.isfile(fastafile) else sys.exit("%s can not be found\n" % fastafile)
	newfasta_handle = open(newffile, "w")
	for record in SeqIO.parse(origfasta, "fasta"):
		if re.search("^\S+;size=\d+;", record.id):
			record.id = re.search("(\S+);size=\d+;", record.id)
			oldrecord = record.id.group(1)
			try:
				if latestsequencesheaders[oldrecord] != None or latestsequencesheaders[oldrecord] != "":
					i = 0
					while i < nrseqheaders[oldrecord]:
						record.id = latestsequencesheaders[oldrecord][i]
						record.description = record.id
						SeqIO.write(record, newfasta_handle, "fasta")
						i += 1
			except KeyError:
				pass
	origfasta.close()
	newfasta_handle.close()

## Degapping FASTA file routine
def DegappingFASTAFile(fastafile):

	if os.path.isfile(fastafile):

		# Naming the output fasta file
		outputfasta = fastafile.replace(".fasta", ".mod.fasta")

		# Starting changing headers to be functional in the pipeline
		origfasta = open(fastafile, "r")
		modfasta = open(outputfasta, "w")
		for record in SeqIO.parse(origfasta, "fasta"):
			record.seq = record.seq.ungap("-")
			SeqIO.write(record, modfasta, "fasta")
		origfasta.close()
		modfasta.close()

		# Replacing the old file for the new one
		os.rename(outputfasta, fastafile)

	else:
		sys.exit("The program cannot open %s\n" % fastafile)

## GNUv3 "window"
def GNUv3():
	sys.exit('''
                    GNU GENERAL PUBLIC LICENSE
                       Version 3, 29 June 2007

 Copyright (C) 2007 Free Software Foundation, Inc. <http://fsf.org/>
 Everyone is permitted to copy and distribute verbatim copies of this license document, but changing it is not allowed.

                            Preamble

  The GNU General Public License is a free, copyleft license for software and other kinds of works.

  The licenses for most software and other practical works are designed to take away your freedom to share and change the works.  By contrast, the GNU General Public License is intended to guarantee your freedom to share and change all versions of a program--to make sure it remains free software for all its users.  We, the Free Software Foundation, use the GNU General Public License for most of our software; it applies also to any other work released this way by its authors.  You can apply it to your programs, too.

  When we speak of free software, we are referring to freedom, not price. Our General Public Licenses are designed to make sure that you have the freedom to distribute copies of free software (and charge for them if you wish), that you receive source code or can get it if you want it, that you can change the software or use pieces of it in new free programs, and that you know you can do these things.

  To protect your rights, we need to prevent others from denying you these rights or asking you to surrender the rights.  Therefore, you have certain responsibilities if you distribute copies of the software, or if you modify it: responsibilities to respect the freedom of others.

  For example, if you distribute copies of such a program, whether gratis or for a fee, you must pass on to the recipients the same freedoms that you received.  You must make sure that they, too, receive or can get the source code.  And you must show them these terms so they know their rights.

  Developers that use the GNU GPL protect your rights with two steps: (1) assert copyright on the software, and (2) offer you this License giving you legal permission to copy, distribute and/or modify it.

  For the developers' and authors' protection, the GPL clearly explains that there is no warranty for this free software.  For both users' and authors' sake, the GPL requires that modified versions be marked as changed, so that their problems will not be attributed erroneously to authors of previous versions.

  Some devices are designed to deny users access to install or run modified versions of the software inside them, although the manufacturer can do so.  This is fundamentally incompatible with the aim of protecting users' freedom to change the software.  The systematic pattern of such abuse occurs in the area of products for individuals to use, which is precisely where it is most unacceptable.  Therefore, we have designed this version of the GPL to prohibit the practice for those products.  If such problems arise substantially in other domains, we stand ready to extend this provision to those domains in future versions of the GPL, as needed to protect the freedom of users.

  Finally, every program is threatened constantly by software patents. States should not allow patents to restrict development and use of software on general-purpose computers, but in those that do, we wish to avoid the special danger that patents applied to a free program could make it effectively proprietary.  To prevent this, the GPL assures that patents cannot be used to render the program non-free.

  The precise terms and conditions for copying, distribution and modification follow.

                       TERMS AND CONDITIONS

  0. Definitions.

  'This License' refers to version 3 of the GNU General Public License.

  'Copyright' also means copyright-like laws that apply to other kinds of works, such as semiconductor masks.

  'The Program' refers to any copyrightable work licensed under this License.  Each licensee is addressed as 'you'. 'Licensees' and 'recipients' may be individuals or organizations.

  To 'modify' a work means to copy from or adapt all or part of the work in a fashion requiring copyright permission, other than the making of an exact copy.  The resulting work is called a 'modified version' of the earlier work or a work 'based on' the earlier work.

  A 'covered work' means either the unmodified Program or a work based on the Program.

  To 'propagate' a work means to do anything with it that, without permission, would make you directly or secondarily liable for infringement under applicable copyright law, except executing it on a computer or modifying a private copy.  Propagation includes copying, distribution (with or without modification), making available to the public, and in some countries other activities as well.

  To 'convey' a work means any kind of propagation that enables other parties to make or receive copies.  Mere interaction with a user through a computer network, with no transfer of a copy, is not conveying.

  An interactive user interface displays 'Appropriate Legal Notices' to the extent that it includes a convenient and prominently visible feature that (1) displays an appropriate copyright notice, and (2) tells the user that there is no warranty for the work (except to the extent that warranties are provided), that licensees may convey the work under this License, and how to view a copy of this License.  If the interface presents a list of user commands or options, such as a menu, a prominent item in the list meets this criterion.

  1. Source Code.

  The 'source code' for a work means the preferred form of the work for making modifications to it. 'Object code' means any non-source form of a work.

  A 'Standard Interface' means an interface that either is an official standard defined by a recognized standards body, or, in the case of interfaces specified for a particular programming language, one that is widely used among developers working in that language.

  The 'System Libraries' of an executable work include anything, other than the work as a whole, that (a) is included in the normal form of packaging a Major Component, but which is not part of that Major Component, and (b) serves only to enable use of the work with that Major Component, or to implement a Standard Interface for which an implementation is available to the public in source code form.  A 'Major Component', in this context, means a major essential component (kernel, window system, and so on) of the specific operating system (if any) on which the executable work runs, or a compiler used to produce the work, or an object code interpreter used to run it.

  The 'Corresponding Source' for a work in object code form means all the source code needed to generate, install, and (for an executable work) run the object code and to modify the work, including scripts to control those activities.  However, it does not include the work's System Libraries, or general-purpose tools or generally available free programs which are used unmodified in performing those activities but which are not part of the work.  For example, Corresponding Source includes interface definition files associated with source files for the work, and the source code for shared libraries and dynamically linked subprograms that the work is specifically designed to require, such as by intimate data communication or control flow between those subprograms and other parts of the work.

  The Corresponding Source need not include anything that users can regenerate automatically from other parts of the Corresponding Source.

  The Corresponding Source for a work in source code form is that same work.

  2. Basic Permissions.

  All rights granted under this License are granted for the term of copyright on the Program, and are irrevocable provided the stated conditions are met.  This License explicitly affirms your unlimited permission to run the unmodified Program.  The output from running a covered work is covered by this License only if the output, given its content, constitutes a covered work.  This License acknowledges your rights of fair use or other equivalent, as provided by copyright law.

  You may make, run and propagate covered works that you do not convey, without conditions so long as your license otherwise remains in force.  You may convey covered works to others for the sole purpose of having them make modifications exclusively for you, or provide you with facilities for running those works, provided that you comply with the terms of this License in conveying all material for which you do not control copyright.  Those thus making or running the covered works for you must do so exclusively on your behalf, under your direction and control, on terms that prohibit them from making any copies of your copyrighted material outside their relationship with you.

  Conveying under any other circumstances is permitted solely under the conditions stated below.  Sublicensing is not allowed; section 10 makes it unnecessary.

  3. Protecting Users' Legal Rights From Anti-Circumvention Law.

  No covered work shall be deemed part of an effective technological measure under any applicable law fulfilling obligations under article 11 of the WIPO copyright treaty adopted on 20 December 1996, or similar laws prohibiting or restricting circumvention of such measures.

  When you convey a covered work, you waive any legal power to forbid circumvention of technological measures to the extent such circumvention is effected by exercising rights under this License with respect to the covered work, and you disclaim any intention to limit operation or modification of the work as a means of enforcing, against the work's users, your or third parties' legal rights to forbid circumvention of technological measures.

  4. Conveying Verbatim Copies.

  You may convey verbatim copies of the Program's source code as you receive it, in any medium, provided that you conspicuously and appropriately publish on each copy an appropriate copyright notice; keep intact all notices stating that this License and any non-permissive terms added in accord with section 7 apply to the code; keep intact all notices of the absence of any warranty; and give all recipients a copy of this License along with the Program.

  You may charge any price or no price for each copy that you convey, and you may offer support or warranty protection for a fee.

  5. Conveying Modified Source Versions.

  You may convey a work based on the Program, or the modifications to produce it from the Program, in the form of source code under the terms of section 4, provided that you also meet all of these conditions:

    a) The work must carry prominent notices stating that you modified it, and giving a relevant date.

    b) The work must carry prominent notices stating that it is released under this License and any conditions added under section 7.  This requirement modifies the requirement in section 4 to 'keep intact all notices'.

    c) You must license the entire work, as a whole, under this License to anyone who comes into possession of a copy.  This License will therefore apply, along with any applicable section 7 additional terms, to the whole of the work, and all its parts, regardless of how they are packaged.  This License gives no permission to license the work in any other way, but it does not invalidate such permission if you have separately received it.

    d) If the work has interactive user interfaces, each must display Appropriate Legal Notices; however, if the Program has interactive interfaces that do not display Appropriate Legal Notices, your work need not make them do so.

  A compilation of a covered work with other separate and independent works, which are not by their nature extensions of the covered work, and which are not combined with it such as to form a larger program, in or on a volume of a storage or distribution medium, is called an 'aggregate' if the compilation and its resulting copyright are not used to limit the access or legal rights of the compilation's users beyond what the individual works permit.  Inclusion of a covered work in an aggregate does not cause this License to apply to the other parts of the aggregate.

  6. Conveying Non-Source Forms.

  You may convey a covered work in object code form under the terms of sections 4 and 5, provided that you also convey the machine-readable Corresponding Source under the terms of this License, in one of these ways:

    a) Convey the object code in, or embodied in, a physical product (including a physical distribution medium), accompanied by the Corresponding Source fixed on a durable physical medium customarily used for software interchange.

    b) Convey the object code in, or embodied in, a physical product (including a physical distribution medium), accompanied by a written offer, valid for at least three years and valid for as long as you offer spare parts or customer support for that product model, to give anyone who possesses the object code either (1) a copy of the Corresponding Source for all the software in the product that is covered by this License, on a durable physical medium customarily used for software interchange, for a price no more than your reasonable cost of physically performing this conveying of source, or (2) access to copy the Corresponding Source from a network server at no charge.

    c) Convey individual copies of the object code with a copy of the written offer to provide the Corresponding Source.  This alternative is allowed only occasionally and noncommercially, and only if you received the object code with such an offer, in accord with subsection 6b.

    d) Convey the object code by offering access from a designated place (gratis or for a charge), and offer equivalent access to the Corresponding Source in the same way through the same place at no further charge.  You need not require recipients to copy the Corresponding Source along with the object code.  If the place to copy the object code is a network server, the Corresponding Source may be on a different server (operated by you or a third party) that supports equivalent copying facilities, provided you maintain clear directions next to the object code saying where to find the Corresponding Source.  Regardless of what server hosts the Corresponding Source, you remain obligated to ensure that it is available for as long as needed to satisfy these requirements.

    e) Convey the object code using peer-to-peer transmission, provided you inform other peers where the object code and Corresponding Source of the work are being offered to the general public at no charge under subsection 6d.

  A separable portion of the object code, whose source code is excluded from the Corresponding Source as a System Library, need not be included in conveying the object code work.

  A 'User Product' is either (1) a 'consumer product', which means any tangible personal property which is normally used for personal, family, or household purposes, or (2) anything designed or sold for incorporation into a dwelling.  In determining whether a product is a consumer product, doubtful cases shall be resolved in favor of coverage.  For a particular product received by a particular user, 'normally used' refers to a typical or common use of that class of product, regardless of the status of the particular user or of the way in which the particular user actually uses, or expects or is expected to use, the product.  A product is a consumer product regardless of whether the product has substantial commercial, industrial or non-consumer uses, unless such uses represent the only significant mode of use of the product.

  'Installation Information' for a User Product means any methods, procedures, authorization keys, or other information required to install and execute modified versions of a covered work in that User Product from a modified version of its Corresponding Source.  The information must suffice to ensure that the continued functioning of the modified object code is in no case prevented or interfered with solely because modification has been made.

  If you convey an object code work under this section in, or with, or specifically for use in, a User Product, and the conveying occurs as part of a transaction in which the right of possession and use of the User Product is transferred to the recipient in perpetuity or for a fixed term (regardless of how the transaction is characterized), the Corresponding Source conveyed under this section must be accompanied by the Installation Information.  But this requirement does not apply if neither you nor any third party retains the ability to install modified object code on the User Product (for example, the work has been installed in ROM).

  The requirement to provide Installation Information does not include a requirement to continue to provide support service, warranty, or updates for a work that has been modified or installed by the recipient, or for the User Product in which it has been modified or installed.  Access to a network may be denied when the modification itself materially and adversely affects the operation of the network or violates the rules and protocols for communication across the network.

  Corresponding Source conveyed, and Installation Information provided, in accord with this section must be in a format that is publicly documented (and with an implementation available to the public in source code form), and must require no special password or key for unpacking, reading or copying.

  7. Additional Terms.

  'Additional permissions' are terms that supplement the terms of this License by making exceptions from one or more of its conditions. Additional permissions that are applicable to the entire Program shall be treated as though they were included in this License, to the extent that they are valid under applicable law.  If additional permissions apply only to part of the Program, that part may be used separately under those permissions, but the entire Program remains governed by this License without regard to the additional permissions.

  When you convey a copy of a covered work, you may at your option remove any additional permissions from that copy, or from any part of it.  (Additional permissions may be written to require their own removal in certain cases when you modify the work.)  You may place additional permissions on material, added by you to a covered work, for which you have or can give appropriate copyright permission.

  Notwithstanding any other provision of this License, for material you add to a covered work, you may (if authorized by the copyright holders of that material) supplement the terms of this License with terms:

    a) Disclaiming warranty or limiting liability differently from the terms of sections 15 and 16 of this License; or

    b) Requiring preservation of specified reasonable legal notices or author attributions in that material or in the Appropriate Legal Notices displayed by works containing it; or

    c) Prohibiting misrepresentation of the origin of that material, or requiring that modified versions of such material be marked in reasonable ways as different from the original version; or

    d) Limiting the use for publicity purposes of names of licensors or authors of the material; or

    e) Declining to grant rights under trademark law for use of some trade names, trademarks, or service marks; or

    f) Requiring indemnification of licensors and authors of that material by anyone who conveys the material (or modified versions of it) with contractual assumptions of liability to the recipient, for any liability that these contractual assumptions directly impose on those licensors and authors.

  All other non-permissive additional terms are considered 'further restrictions' within the meaning of section 10.  If the Program as you received it, or any part of it, contains a notice stating that it is governed by this License along with a term that is a further restriction, you may remove that term.  If a license document contains a further restriction but permits relicensing or conveying under this License, you may add to a covered work material governed by the terms of that license document, provided that the further restriction does not survive such relicensing or conveying.

  If you add terms to a covered work in accord with this section, you must place, in the relevant source files, a statement of the additional terms that apply to those files, or a notice indicating where to find the applicable terms.

  Additional terms, permissive or non-permissive, may be stated in the form of a separately written license, or stated as exceptions; the above requirements apply either way.

  8. Termination.

  You may not propagate or modify a covered work except as expressly provided under this License.  Any attempt otherwise to propagate or modify it is void, and will automatically terminate your rights under this License (including any patent licenses granted under the third paragraph of section 11).

  However, if you cease all violation of this License, then your license from a particular copyright holder is reinstated (a) provisionally, unless and until the copyright holder explicitly and finally terminates your license, and (b) permanently, if the copyright holder fails to notify you of the violation by some reasonable means prior to 60 days after the cessation.

  Moreover, your license from a particular copyright holder is reinstated permanently if the copyright holder notifies you of the violation by some reasonable means, this is the first time you have received notice of violation of this License (for any work) from that copyright holder, and you cure the violation prior to 30 days after your receipt of the notice.

  Termination of your rights under this section does not terminate the licenses of parties who have received copies or rights from you under this License.  If your rights have been terminated and not permanently reinstated, you do not qualify to receive new licenses for the same material under section 10.

  9. Acceptance Not Required for Having Copies.

  You are not required to accept this License in order to receive or run a copy of the Program.  Ancillary propagation of a covered work occurring solely as a consequence of using peer-to-peer transmission to receive a copy likewise does not require acceptance.  However, nothing other than this License grants you permission to propagate or modify any covered work.  These actions infringe copyright if you do not accept this License.  Therefore, by modifying or propagating a covered work, you indicate your acceptance of this License to do so.

  10. Automatic Licensing of Downstream Recipients.

  Each time you convey a covered work, the recipient automatically receives a license from the original licensors, to run, modify and propagate that work, subject to this License.  You are not responsible for enforcing compliance by third parties with this License.

  An 'entity transaction' is a transaction transferring control of an organization, or substantially all assets of one, or subdividing an organization, or merging organizations.  If propagation of a covered work results from an entity transaction, each party to that transaction who receives a copy of the work also receives whatever licenses to the work the party's predecessor in interest had or could give under the previous paragraph, plus a right to possession of the Corresponding Source of the work from the predecessor in interest, if the predecessor has it or can get it with reasonable efforts.

  You may not impose any further restrictions on the exercise of the rights granted or affirmed under this License.  For example, you may not impose a license fee, royalty, or other charge for exercise of rights granted under this License, and you may not initiate litigation (including a cross-claim or counterclaim in a lawsuit) alleging that any patent claim is infringed by making, using, selling, offering for sale, or importing the Program or any portion of it.

  11. Patents.

  A 'contributor' is a copyright holder who authorizes use under this License of the Program or a work on which the Program is based.  The work thus licensed is called the contributor's 'contributor version'.

  A contributor's 'essential patent claims' are all patent claims owned or controlled by the contributor, whether already acquired or hereafter acquired, that would be infringed by some manner, permitted by this License, of making, using, or selling its contributor version, but do not include claims that would be infringed only as a consequence of further modification of the contributor version.  For purposes of this definition, 'control' includes the right to grant patent sublicenses in a manner consistent with the requirements of this License.

  Each contributor grants you a non-exclusive, worldwide, royalty-free patent license under the contributor's essential patent claims, to make, use, sell, offer for sale, import and otherwise run, modify and propagate the contents of its contributor version.

  In the following three paragraphs, a 'patent license' is any express agreement or commitment, however denominated, not to enforce a patent (such as an express permission to practice a patent or covenant not to sue for patent infringement).  To 'grant' such a patent license to a party means to make such an agreement or commitment not to enforce a patent against the party.

  If you convey a covered work, knowingly relying on a patent license, and the Corresponding Source of the work is not available for anyone to copy, free of charge and under the terms of this License, through a publicly available network server or other readily accessible means, then you must either (1) cause the Corresponding Source to be so available, or (2) arrange to deprive yourself of the benefit of the patent license for this particular work, or (3) arrange, in a manner consistent with the requirements of this License, to extend the patent license to downstream recipients. 'Knowingly relying' means you have actual knowledge that, but for the patent license, your conveying the covered work in a country, or your recipient's use of the covered work in a country, would infringe one or more identifiable patents in that country that you have reason to believe are valid.

  If, pursuant to or in connection with a single transaction or arrangement, you convey, or propagate by procuring conveyance of, a covered work, and grant a patent license to some of the parties receiving the covered work authorizing them to use, propagate, modify or convey a specific copy of the covered work, then the patent license you grant is automatically extended to all recipients of the covered work and works based on it.

  A patent license is 'discriminatory' if it does not include within the scope of its coverage, prohibits the exercise of, or is conditioned on the non-exercise of one or more of the rights that are specifically granted under this License.  You may not convey a covered work if you are a party to an arrangement with a third party that is in the business of distributing software, under which you make payment to the third party based on the extent of your activity of conveying the work, and under which the third party grants, to any of the parties who would receive the covered work from you, a discriminatory patent license (a) in connection with copies of the covered work conveyed by you (or copies made from those copies), or (b) primarily for and in connection with specific products or compilations that contain the covered work, unless you entered into that arrangement, or that patent license was granted, prior to 28 March 2007.

  Nothing in this License shall be construed as excluding or limiting any implied license or other defenses to infringement that may otherwise be available to you under applicable patent law.

  12. No Surrender of Others' Freedom.

  If conditions are imposed on you (whether by court order, agreement or otherwise) that contradict the conditions of this License, they do not excuse you from the conditions of this License.  If you cannot convey a covered work so as to satisfy simultaneously your obligations under this License and any other pertinent obligations, then as a consequence you may not convey it at all.  For example, if you agree to terms that obligate you to collect a royalty for further conveying from those to whom you convey the Program, the only way you could satisfy both those terms and this License would be to refrain entirely from conveying the Program.

  13. Use with the GNU Affero General Public License.

  Notwithstanding any other provision of this License, you have permission to link or combine any covered work with a work licensed under version 3 of the GNU Affero General Public License into a single combined work, and to convey the resulting work.  The terms of this License will continue to apply to the part which is the covered work, but the special requirements of the GNU Affero General Public License, section 13, concerning interaction through a network will apply to the combination as such.

  14. Revised Versions of this License.

  The Free Software Foundation may publish revised and/or new versions of the GNU General Public License from time to time.  Such new versions will be similar in spirit to the present version, but may differ in detail to address new problems or concerns.

  Each version is given a distinguishing version number.  If the Program specifies that a certain numbered version of the GNU General Public License 'or any later version' applies to it, you have the option of following the terms and conditions either of that numbered version or of any later version published by the Free Software Foundation.  If the Program does not specify a version number of the GNU General Public License, you may choose any version ever published by the Free Software Foundation.

  If the Program specifies that a proxy can decide which future versions of the GNU General Public License can be used, that proxy's public statement of acceptance of a version permanently authorizes you to choose that version for the Program.

  Later license versions may give you additional or different permissions.  However, no additional obligations are imposed on any author or copyright holder as a result of your choosing to follow a later version.

  15. Disclaimer of Warranty.

  THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW.  EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM 'AS IS' WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.

  16. Limitation of Liability.

  IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MODIFIES AND/OR CONVEYS THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.

  17. Interpretation of Sections 15 and 16.

  If the disclaimer of warranty and limitation of liability provided above cannot be given local legal effect according to their terms, reviewing courts shall apply local law that most closely approximates an absolute waiver of all civil liability in connection with the Program, unless a warranty or assumption of liability accompanies a copy of the Program in return for a fee.

                     END OF TERMS AND CONDITIONS

            How to Apply These Terms to Your New Programs

  If you develop a new program, and you want it to be of the greatest possible use to the public, the best way to achieve this is to make it free software which everyone can redistribute and change under these terms.

  To do so, attach the following notices to the program.  It is safest to attach them to the start of each source file to most effectively state the exclusion of warranty; and each file should have at least the 'copyright' line and a pointer to where the full notice is found.

    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) <year>  <name of author>

    This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

Also add information on how to contact you by electronic and paper mail.

  If the program does terminal interaction, make it output a short notice like this when it starts in an interactive mode:

    <program>  Copyright (C) <year>  <name of author>
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.

The hypothetical commands `show w' and `show c' should show the appropriate parts of the General Public License.  Of course, your program's commands might be different; for a GUI interface, you would use an 'about box'.

  You should also get your employer (if you work as a programmer) or school, if any, to sign a 'copyright disclaimer' for the program, if necessary. For more information on this, and how to apply and follow the GNU GPL, see <http://www.gnu.org/licenses/>.

  The GNU General Public License does not permit incorporating your program into proprietary programs.  If your program is a subroutine library, you may consider it more useful to permit linking proprietary applications with the library.  If this is what you want to do, use the GNU Lesser General Public License instead of this License.  But first, please read <http://www.gnu.org/philosophy/why-not-lgpl.html>.
\n''')

## Help "window"
def HelpingU():
	sys.exit('''
NAME: Quantification of Representative Sequence (QRS.py)

AUTHOR: Enrique Gonzalez-Tortuero <gonzalez@igb-berlin.de>

DESCRIPTION: QRS is a script (originally written in Perl 5.14.2 but now in Python 3.4) that allows analysing NGS datasets automatically to study population structure. Currently, this pipeline is optimized for pyrosequencing platform and Illumina sequences but it is possible to use it for another NGS technologies like Ion Torrent. This program may operate either in batch processing mode where all sequences are automatically analysed in an unsupervised way or it may interact with the user at various checkpoints if no parameters have been specified prior to execution.

	USAGE:

The program QRS can be executed as batch command or as an interactive program if no parameters are defined. To understand how QRS works in batch mode, we provide several use-case scenarios.

	Creating a reference HMM file (-RM)
	-----------------------------------

# Basic parameters

When this option is used, the batch command can be executed using only basic parameters:
		QRS.py -RM -folder=test -reffiles=test.fasta

In this case, we call QRS to create a reference HMM file of test.fasta in the folder test using default parameters (the aligner is PRANK and QRS will call ReformAl to post-process the alignment).
However, QRS can create more than a single HMM file if you specify different reference files joined with a single dash (-) in the -reffiles parameter, like in the following example, where the program will create two HMM files using default parameters:
		QRS.py -RM -folder=test -reffiles=test1.fasta-test2.fasta

# Alignment options

Another option that you can modify is the multiple sequence alignment program. You can choose the program with the -aligner parameter. PRANK is the default aligner program, like in this example where QRS will create a single reference HMM file:
		QRS.py -RM -folder=test -aligner=prank -reffiles=test.fasta

To employ a different aligner for the sequence alignment task simply modify the -aligner parameter. In the following example, QRS will call MUSCLE to align a reference file to create a single reference HMM file:
		QRS.py -RM -folder=test -aligner=muscle -reffiles=test.fasta

Here, we use MAFFT with an iterative global alignment algorithm called G-INS-i to perform an accurate alignment:
		QRS.py -RM -folder=test -aligner=mafft-ginsi -reffiles=test.fasta

The last option you may modify is whether to use ReformAl to post-process the alignment or not. By default, this option is enabled (-reformal=yes), like in the following example, where QRS will create two reference HMM files after calling GramAlign to align your reference sequences:
		QRS.py -RM -folder=test -aligner=gramalign -reformal=yes -reffiles=test1.fasta-test2.fasta

If you want to disable this option, you have to type -reformal=no. In this example, QRS will align your three reference files using the default aligner (i.e., PRANK) and will not curate the alignments before creating three reference HMM files:
		QRS.py -RM -folder=test -reformal=no -reffiles=test1.fasta-test2.fasta-test3.fasta

In the following example, we call QRS to align four reference files using MAFFT with their iterative local alignment algorithm (L-INS-i) and we do not want to post-process the alignments:
		QRS.py -RM -folder=test -aligner=mafft-linsi -reformal=no -reffiles=test1.fasta-test2.fasta-test3.fasta-test4.fasta

Finally, if you are working with standard alignments (like SILVA [REFERENCE] or GreenGenes [REFERENCE]), perhaps you want to avoid the alignment step. To do it, you can use the parameter "-noalign" like in the following example, where we want to create a HMM file from SILVA alignment (in fasta format):
		QRS.py -RM -folder=test -noalign -reffiles=silva.fasta

	Describing the NGS data set (-AM1)
	----------------------------------

# Basic parameters

When this option is used, the batch command cannot be executed using only basic parameters. Instead, you have to specify if you have a NGS FASTA file or a NGS FASTQ file. In this example the input file is a FASTA file: 
	QRS.py -AM1 -folder=test -informat=fasta -fasta=file.fna

Here, the input file is a FASTQ file:
	QRS.py -AM1 -folder=test -informat=fastq -fastq=file.fastq

# Another input files
If you have a FASTA file as input, you can add a quality file to do all basic statistics on the base quality using the parameter -quality:
	QRS.py -AM1 -folder=test -informat=fasta -fasta=file.fna -quality=file.qual

If you have a paired FASTQ file, you have to add the parameter -paired=yes and the name of the paired FASTQ:
	QRS.py -AM1 -folder=test -informat=fastq -fastq=file.fastq -paired=yes -fastq2=file2.fastq

	Processing a NGS data set (-AM2)
	--------------------------------

# Basic parameters

Before executing this part, you have to evaluate the characteristics of your data set according to the previous step (see AM1 section). As all NGS datasets are different depending on the case study, there are no default parameters to filter and trim sequences according to length, quality and complexity. 
However, in a hypothetical case that you do not need to filter your data set, the input parameters are the input files added as describe before (see AM1 section). Moreover, you have to add the reference HMM file (created in -RM step, see RM section for more details), oligos and design files (read below to know more about these files) and specify if you have reverse barcoded primers or not. In the following example, we use a NGS FASTA file with quality file and our data set has reverse barcoded primers (-bry): 
	QRS.py -AM2 -folder=test -informat=fasta -fasta=file.fna -quality=file.qual -hmmfile=test.hmm -oligos=oligos.csv -design=design.csv -bry

The oligos file is a plain text file that has two columns splitted by tabs (except if the label is barcode, where there are a third column). The first column contains always a label that indicates forward or reverse adapter, forward or reverse primer and barcode. The second contains the DNA sequence for each element and, in the case of barcode label, the third one indicates the barcode ID (MID1, MID2, MID3...). An example of this file is as follows:
		   seqadapfor	SEQUENCE
		   seqadaprev	SEQUENCE
		   forward	SEQUENCE
		   reverse	SEQUENCE
		   barcode	SEQUENCE	BARID

The design file is also a plain text file that has two or three columns (depending of the use of reverse barcoded primers) splitted by tabs. An example of these file is as follows:
		   BARID1	(BARID2)	SAMPLEID

Here, BARID1 is the Barcode ID for the forward primer, BARID2 is the Barcode ID for the reverse primer (if exists) and SAMPLEID is the name of the sample.

# Filtering by length

As it is said in AM2.1 section, there are no default parameters to filter and trim sequences according to length, quality and complexity because it depends of your data set. In the following example calls QRS to accept sequences that are greater than 300 bp in a NGS FASTQ file and have no reverse barcoded primers:
	QRS.py -AM2 -folder=test -informat=fastq -fastq=file.fastq -hmmfile=test.hmm -minlen=300 -oligos=oligos.csv -design=design.csv -brn

In this example, QRS is executed to accept sequences that have between 355 and 500 bp (the data is given as NGS FASTA file without quality file):
	QRS.py -AM2 -folder=test -informat=fasta -fasta=file.fna -hmmfile=test.hmm -minlen=355 -maxlen=500 -oligos=oligos.csv -design=design.csv -bry

# Filtering by GC content

Another possible filter is based on the GC content. In the following example QRS is called to accept sequences that have more than 60% GC content in a NGS FASTQ file and have no reverse barcoded primers:
	QRS.py -AM2 -folder=test -informat=fastq -fastq=file.fastq -hmmfile=test.hmm -mingc=60 -oligos=oligos.csv -design=design.csv -brn

In this example, QRS is executed to accept sequences that have between 40 and 50% GC content (the data is given as NGS FASTA file without quality file):
	QRS.py -AM2 -folder=test -informat=fasta -fasta=file.fna -hmmfile=test.hmm -mingc=40 -maxgc=50 -oligos=oligos.csv -design=design.csv -bry

# Filtering and trimming according to quality

We offer another filter based on base quality. In this example, QRS accepts only sequences that have at least a mean sequence quality value of 25:
	QRS.py -AM2 -folder=test -informat=fastq -fastq=file.fastq -hmmfile=test.hmm -minqual=25 -oligos=oligos.csv -design=design.csv -bry

If you have a FASTA file, it is convenient to have a quality file to filter by base quality. In the following example, QRS will filter an input NGS FASTA file by length (accepting only sequences that are between 375 and 480 bp long) and quality (rejecting sequences that have a mean sequence quality score smaller than 28). This data set has no reverse barcoded primers:
	QRS.py -AM2 -folder=test -informat=fasta -fasta=file.fna -quality=file.qual -hmmfile=test.hmm -oligos=oligos.csv -design=design.csv -minlen=375 -maxlen=480 -minqual=28 -brn

The following option is used to trim bases that have an insufficient quality score. In the following example, QRS will trim all nucleotides with a quality value less than 25 in the beginning and at the end of the sequences in a non-reverse barcoded NGS data set:
	QRS.py -AM2 -folder=test -informat=fasta -fasta=file.fna -quality=file.qual -hmmfile=test.hmm -oligos=oligos.csv -design=design.csv -trimqualleft=25 -trimqualright=25 -brn

Finally, another way to filter your dataset according to quality is allowing ambiguous characters (Ns). The number of allowed Ns in your dataset can be modified with the parameter -allowns. In the following example, QRS accepts all sequences that have at most one ambiguous character in the sequence:
	QRS.py -AM2 -folder=folder -informat=fastq -fastq=file.fastq -hmmfile=ref.hmm -oligos=oligos.csv -design=design.csv -allowns=1 -bry

However, it is not recommended to allow ambiguous characters as they might indicate bad quality nucleotides (Huse et al. 2007). In the following example, QRS removes all sequences that have more than one ambiguous character:
	QRS.py -AM2 -folder=folder -informat=fasta -fasta=file.fna -hmmfile=ref.hmm -oligos=oligos.csv -design=design.csv -allowns=0 -bry

# Filtering by complexity

Finally, if you want to remove all low complexity sequences like homopolymers, you can do this with-filmet and -filthr parameters. The first argument defines the method to remove this kind of sequences (DUST (Morgulis et al. 2006) or Entropy-based filter) and the second argument defines the threshold for the filtering method. For more details on these filtering methods, see PrinSeq manual (http://prinseq.sourceforge.net/manual.html#QCCOMPLEXITY). In the following example, we use DUST to remove all sequences with a complexity greater than 7:
	QRS.py -AM2 -folder=folder -informat=fasta -fasta=file.fna -hmmfile=ref.hmm -oligos=oligos.csv -design=design.csv -filmet=dust -filthr=7 -bry

In this example, we use entropy to remove all sequences with a complexity smaller than 70:
	QRS.py -AM2 -folder=folder -informat=FASTA -FASTA=file.fna -hmmfile=ref.hmm -oligos=oligos.csv -design=design.csv -filmet=entropy -filthr=70 -bry

In this example, a NGS dataset (that consists in a FASTQ file) is filtered by length (accepting only sequences that are between 390 and 500 bp long), base quality (removing all sequences with mean sequence quality score less than 28), and low-complexity based on the DUST filter (considering all sequences with values greater than 5 as low complexity sequences):
	QRS.py -AM2 -folder=folder -informat=fasta -fasta=file.fna -quality=file.qual -hmmfile=ref.hmm -oligos=oligos.csv -design=design.csv -minlen=390 -maxlen=500 -minqual=28 -filmet=dust -filthr=5 -bry

Finally, it is feasible to trim poly-A/Ts tails. In this example, QRS trims the regions that have more than three followed A/T at the beginning or the end of the sequences:
	QRS.py -AM2 -folder=test -informat=fastq -FASTQ=file.fastq -hmmfile=test.hmm -trimtails=3 -oligos=oligos.csv -design=design.csv -bry

# Classifying sequences as a specific marker

You can also modify the maximum allowed e-value to classify a sequence as a specific marker using a HMM profile (-hmmthr). By default, this argument is set to 10-10 as then similar results as with the usage of BLAST are achieved (Altschul et al. 1990). We do not recommend changing this value. However, you can change the value like in the following example where QRS classifies a sequence as a specific marker if the e-value is smaller than 10-25:
	QRS.py -AM2 -folder=folder -informat=fasta -fasta=file.fna -hmmfile=ref.hmm -oligos=oligos.csv -design=design.csv -hmmthr=1E-25 -bry

# Assigning sequences to samples

If you want to modify the maximum number of allowed mismatches to detect the barcodes in your data set, you can use the parameter -allowmis. In this example, QRS classifies a non reverse barcoded data set in different samples with a maximum number of allowed mismatches from 1 to positive infinite. By default, the value is 1 and we do not recommend to modify this value, although it is possible to change it like in this example, where QRS considers 2 mismatches to detect the barcodes in a reverse barcoded data set:
	QRS.py -AM2 -folder=folder -informat=FASTQ -FASTQ=file.FASTQ -hmmfile=ref.hmm -oligos=oligos.csv -design=design.csv -allowmis=2 -bry

# Clustering step

Another parameter you can modify is -cutoff. This parameter is used to cluster sequences in order to de-noise your data set, i. e., to remove all spurious nucleotides across all sequences. By default, QRS calculates this value according to the CD-HIT-OTU algorithm (Li et al. 2012) but you can define the value by a number between 0.00 and 1.00. For example, if you want to cluster all sequences that have 99.5% of similarity, you can type:
	QRS.py -AM2 -folder=folder -informat=fasta -fasta=file.fna -hmmfile=ref.hmm -oligos=oligos.csv -design=design.csv -cutoff=0.995 -brn

The parameter -minclustersize specifies the minimum size of cluster to consider that analysed sequences are not sequencing artifacts. By default, this filter is enabled and all clusters that have at least three sequences are considered as good clusters. If you want to disable this argument, you have to type -minclustersize=0, like in this example:
	QRS.py -AM2 -folder=folder -informat=fasta -fasta=file.fna -hmmfile=ref.hmm -oligos=oligos.csv -design=design.csv -minclustersize=0 -bry

# Alignment step

The last options you may modify concern the use of the alignment program and ReformAl to post-process the alignment. For more information, see the RM examples.

# A complex example

In the following example, QRS retrieves all sequences that have more than 350 bp, a mean sequence quality score of 30, no ambiguous characters, and a sequence complexity greater than 79 according to the entropy-based filter. This data set has reverse paired barcodes and the script discards all clusters that have only one sequence. Finally, QRS uses KAlign to align all accepted sequences and this alignment is not fine tuned:
	QRS.py -AM2 -folder=folder -informat=fasta -fasta=file.fna -quality=file.qual -hmmfile=ref.hmm -oligos=oligos.csv -design=design.csv -minlen=350 -minqual=30 -allowns=0 -filmet=entropy -filthr=79 -bry -minclusterzise=2 -aligner=kalign -reformal=no

	Classifying and quantifying all NGS sequences into representative sequence variants (-AM3)
	------------------------------------------------------------------------------------------

# Basic parameters

When this option is used, the batch command can be executed with only basic parameters, like in the following example:
	QRS.py -AM3 -folder=test -fasta=myaligneddata.fasta -derep=derep.list -names=names.list -sample=sample.csv -outfile=allsamples

In the previous example, QRS will use myaligneddata alignment to run TCS for the data and then retrieve a TCS-types FASTA file called allsamples.fasta and a frequencies matrix file called allsamples.abs.csv. The program makes use of some information from the samples.csv file. This file which is a plain text file that contains the following information always splitted by tabs:

                SAMPLEID	DATA1	(...)	DATAN

Here, SAMPLEID is the name of the sample and DATA are data you want to put here (Place, year, species...).

derep.list and names.list are files produced in the -AM2 step during the clustering step.

# Statistical parsimony (SP) parameters

In order to use statistical parsimony to pool sequences, you must to use the '-method=SP' parameter:
	QRS.py -AM3 -folder=test -fasta=myaligneddata.fasta -method=SP -derep=derep.list -names=names.list -sample=sample.csv -limitid=0.95 -outfile=allsamples

You can modify the limit threshold that you want to use to determine representative sequences variants with the -limitid parameter. By default, QRS considers that two sequences belong to the same TCS-type if they have more similarity than 0.95, like in the following example:
	QRS.py -AM3 -folder=test -fasta=myaligneddata.fasta -method=SP -derep=derep.list -names=names.list -sample=sample.csv -limitid=0.95 -outfile=allsamples

In this example, the limit threshold to determine if two sequences belong from the same representative sequence variant is 0.995:
	QRS.py -AM3 -folder=test -fasta=myaligneddata.fasta -method=SP -derep=derep.list -names=names.list -sample=sample.csv -limitid=0.995 -outfile=allsamples

Another parameter you can modify is the fact that you want to consider gaps as another character or not (-considergaps). By default, it is enabled (-considergaps=1):
	QRS.py -AM3 -folder=test -fasta=myaligneddata.fasta -method=SP -derep=derep.list -names=names.list -sample=sample.csv -limitid=0.95 -considergaps=1 -outfile=allsamples

However, if you don't want to consider gaps, disable such parameter using -considergaps=0:
	QRS.py -AM3 -folder=test -fasta=myaligneddata.fasta -method=SP -derep=derep.list -names=names.list -sample=sample.csv -limitid=0.95 -considergaps=0 -outfile=allsamples

Finally, the last parameter you want to modify about the Statistical parsimony algorithm is related with the transition/transversion relation (-transition=). As Templeton et al (1992) indicated in their paper, this parameter has only three possible values: 1 (extreme bias to transitions in the DNA - Kimura-2-Parameter model), 2 or 3 (no differences between transitions and transversion - Jukes-Cantor model). By default, this parameter is set as 1:
	QRS.py -AM3 -folder=test -fasta=myaligneddata.fasta -method=SP -sample=sample.csv -limitid=0.95 -considergaps=1 -transition=1 -outfile=allsamples

Although we don't recommend to modify this value, if you consider than the probability of transition is the same than probability of transversion (i.e., Jukes-Cantor model), you should modify this value as follows:
	QRS.py -AM3 -folder=test -fasta=myaligneddata.fasta -method=SP -sample=sample.csv -limitid=0.95 -considergaps=1 -transition=3 -outfile=allsamples

# Neighbour-joining (NJ) parameters

From QRS v. 1.2.0, you can pool your sequences using NJ clustering. To do this, you must to type '-method=NJ' in batch mode. The rest of the parameters are the same as in SP, except that there is no '-considergaps' and '-transition' parameters:
	QRS.py -AM3 -folder=test -fasta=myaligneddata.fasta -method=NJ -sample=sample.csv -derep=derep.list -names=names.list -limitid=0.95 -outfile=allsamples

REQUIREMENTS:

Before using this pipeline, the following Python modules, R libraries and programs should be installed:

* Python modules:
	- Biopython (Bio module; Cock et al. 2009)
	- Numerical Python (numpy module; Oliphant 2007)
	- Processes and System information Utility (psutil module)
	- R from Python (rpy2 module)
	- Scientific Python (scipy module; Oliphant 2007)

* R packages:
	- ape (Paradis et al. 2004)
	- seqinr (Charif & Lobry 2007)

* Programs:
	- PrinSeq (Schmieder & Edwards 2011): it is used to generate summary statistics of sequence and quality data and to filter, reformat and trim NGS data. This program is publicly available at http://prinseq.sourceforge.net under the GPLv3 licence.
	- PEAR (Zhang et al. 2014): it is used to merge all paired fastq files after filtering all your sequences according to length, GC content, quality and complexity with PrinSeq. This program is publicly available at http://www.exelixis-lab.org/web/software/pear under the GPLv3 licence.
	- CutAdapt (Martin 2011): a Python script that removes primers and adapters. This program is publicly available at http://code.google.com/p/cutadapt under the MIT licence.
	- HMMER 3.1b (Finn et al 2011): it is used to search sequences using probabilistic models called profile hidden Markov models (profile HMM). This program is publicly available at http://hmmer.janelia.org/ under GPLv3 licence.
	- USEARCH v. 7.0.1003+ (Edgar 2010): it is used to cluster sequences and detect chimaeras according to the UCHIME algorithm (Edgar et al. 2011). This program is available after requiring it according to the authors' page (http://www.drive5.com/usearch/).
	- At least one of the following aligners:
		* Clustal Omega (Sievers et al. 2011): a hybrid aligner that mixes progressive multiple sequence alignments with the use of Markov models. This program is publicly available at http://www.clustal.org/omega/ under GPLv2 licence.
		* FSA (Bradley et al. 2009): a probabilistic aligner that uses the sequence annealing technique (Schwartz & Pachter 2007) for constructing a multiple alignment from pairwise homology estimation. FSA is publicly available at http://fsa.sourceforge.net/ under GPL licence.
		* GramAlign (Russell et al. 2008): a time-efficient progressive aligner which estimates distances according to the natural grammar present in nucleotide sequences. This program is freely available at http://bioinfo.unl.edu/gramalign.php
		* KAlign (Lassmann & Sonnhammer 2005): a progressive aligner based on Wu-Manber string-matching algorithm. KAlign is publicly available at http://msa.sbc.su.se under GPLv2 licence.
		* MAFFT (Katoh & Standley 2013): an iterative/progressive alignment program that is publicly available at http://mafft.cbrc.jp/alignment/software/ under BSD licence.
		* MUSCLE (Edgar 2004): an iterative aligner that it is freely available at http://www.drive5.com/muscle/.
		* Opal (Wheeler & Kececioglu 2007): a progressive aligner that it is freely available at http://opal.cs.arizona.edu/.
		* PicXAA (Sahraeian & Yoon 2010): a probabilistic non-progressive alignment algorithm that finds multiple sequence alignments with maximum expected accuracy that it is publicly available at http://gsp.tamu.edu/picxaa under GPLv3 licence.
		* PRANK (Loytynoja & Goldman 2005): a probabilistic multiple alignment program based on maximum likelihood methods used in phylogenetics that considers evolutionary distances between sequences. This program is publicly available at http://code.google.com/p/prank-msa under GPLv3 licence.

Although the aforementioned aligners have already been tested with the QRS, other aligners may be added to use in this pipeline. If you want to work with another alignment program, please feel free to contact with the author (qrspipeline@gmail.com) to include it in the source code of QRS. 
Finally, optional recommended requirements might be the installation of ReformAlign (Lyras and Metzler 2014). ReformAlign is a recently proposed profile-based meta-alignment approach that aims to post-process existing alignments via the employment of standard profiles. This program is publicly available at http://evol.bio.lmu.de/_statgen/software/reformalign under GPLv3 licence.

HISTORY: 

v 1.4.0 (DEVEL) - Major issues:
	  * Added the parameter "-noalign" in "-AM2" step. It is only useful if you want to launch the "-AM3" step using NJ clustering.
	  * Added the new parameter "-namesample" that it is mandatory when you have FASTA/FASTQ files for different samples. In this case, this name is used in order to define sequences in each sample and promoting a correct analysis in the "Inferring Representative Sequences" step. The name of the sample should be as is in the first column of the SAMPLES file.
	  * Changed the way of indicating the number of allowed mismatches in order to detect the barcoded tags into the sequences. Now it is referred as the number of base pairs that could be indicate a mismatch in the barcode detection instead of a percentage. Additionally, this step is optional only if you want to assign sequences to samples.
          * When the obtained dataset in "-RM" and "-AM2" steps has more than a thousand sequences, QRS can modify the behaviour of the aligners in order to deal properly with this dataset according to the aligners manuals.
	  * Added optional MAFFT parameters in order to deal with huge datasets: "-aligner=mafft-dpparttree", "-aligner=mafft-fastaparttree" and "-aligner=mafft-parttree". All of them are described in Katoh & Toh (2007) and they are recommended for more than 10,000 sequences according to the authors.
	  * Added a new optional "aligner" that it is the most accurate way to align huge datasets at this moment ("-aligner=muscle-profile"). This method divides the dataset in two: a small subset of ten sequences that belongs to the most abundant sequences ("core") and the rest of the sequences ("rest"). The core subset is aligned using muscle with default parameters and, then, every sequence in the rest subset is aligned against the core using a profile alignment. This method is well recommended when you have several thousands of sequences due to the accuracy of this method (Sievers et al. 2013). However, this method is really slow and it can take more than a day.
	  * Simplified code in order to detect chimeras and de-noising sequences according to USEARCH v. 7.0.1003+ manual.
	  * Added new option in the parameter "-filtermethod=" (in "-AM2" step) that allows to execute HECTOR (Wirawan et al. 2014) in order to trim homopolymers in 454 reads.

	  Minor issues:
	  * Added Opal (Wheeler & Kececioglu 2007) as a new aligner ("-aligner=opal").
	  * In the assigning sequences to samples, you may know in this step hoy many sequences are in each sample after removing primers.
	  * Modified CLUSTAL Omega parameters in order to deal properly fastafiles in "-RM" and "-AM2" step.
	  * Fixed small bug about IDs that are no present in DEREP in NAMES files but are in the fasta file when "-AM3" step is executed.
	  * Fixed small bug about IDs after clustering sequences when a cluster is marked as a chimera.

v 1.3.0 - Major issues:
	  * Fixed a bug in calculation of absolute frequencies using NJ. Now, the log file shows the total number sequences for each network while the absolute matrix should show the number of sequences for each network in every sample.
	  * Fixed a bug in the number of allowed mismatches to detect primers and number of allowed ambiguities (Ns). Now, the parameters "-allowedmis" and "-allowns" should be work properly.
	  * Fixed an error about displaying the help in the program using batch mode. Now, the option "--help" should be work properly.
	  * Fixed an error in Statistical Parsimony method on the "non-parsimony" probability correction. As this step could underestimate the number of networks, we employed the proposed solution of Templeton et al. (1992) to collapse all non-corrected networks only in a single step because it is compute-expensive (and not in all potential steps as they originally proposed). This fix generates the most similar results than TCS (Clement et al. 2000)

v 1.2.0 - Major issues:
          * Added a new parameter "-method" in order to pool your sequences according to Neighbour Joining ("-method=NJ") or Statistical Parsimony ("-method=SP").
	  * Fixed a bug in calculation of absolute frequencies matrix using SP in the last version. Now, the log file shows the total number sequences for each network while the absolute matrix should show the number of sequences for each network in every sample.
	  * Fixed an error in Statistical Parsimony method on the "non-parsimony" probability correction. This step is made after determining the most representative sequences for each network considering only parsimony probability.

	  Minor issues:
          * Added in the output of interactive batch mode the equivalence to the command line options.
	  * Added the verbose output in ReformAlign
	  * Fixed bugs in command line options in -RM options.
	  * Fixed bugs in number of CPUs to run HMMER in -RM options.
	  * Updated references about ReformAlign

v 1.1.0 - Major issues:
          * Script translated to Python 3 with access to R (thanks to RPy2 project)
          * Script infers representative sequence variants according to Statistical Parsimony (Templeton et al. 1992) without executing TCS. This important modification generates the following changes:
                      * Removed "Declustering" step at the end of -AM2 step: In the previous version, the QRS pipeline declustered all sequences after the alignment step because TCS requires a non-clusterized dataset in order to run Statistical Parsimony. As we implemented such algorithm in our program without executing TCS, this step is unnecessary and the output from this step is, directly, a clustered alignment file in order to apply such algorithm. For this reason, the FASTA header must have the following structure: ">(Sample_Name)_(NumberID);size=(Total number of sequences that are identical)".
                      * Added "-limitid", "-considergaps" and "-transition" parameters. All of them are parameter related to Statistical Parsimony according to the Templeton et al. (1992) algorithm. For more info about how to use them, read all -AM3 step examples.
                      * Requiral of DEREP and NAMES files in -AM3 step: such files are produced after clustering sequences in USEARCH and created in QRS pipeline as plain text files as "derep.list" and "names.list". These files contains two columns separated by a tabulation. The first column have the name for the sequence and the second column have all sequences that are considered the same (in the case of derep.list are 100%% identical, while in names.list are synonims due to the denoising step). In batch mode, the parameters are "-derep=(Namefile)" and "-names=(Namefile)", respectively.
          * Added a new process in the -AM2 step between the filtering sequences by length, quality, GC content and complexity using PrinSeq and the filtering sequences by a specific marker using HMMER processes. This new process is called "Merging sequences using PEAR" and it is only useful if you are running paired fastq files. PEAR is a required program in order to run this pipeline (see REQUIREMENTS).

	  Minor issues:
          * Added "-noalign" parameter in the -RM step. It is only useful if you are working with standard alignments like SILVA or GreenGenes.
          * Added a new aligner PicXAA (Sahraeian & Yoon 2010). In order to use it, you can choose two options: "-aligner=picxaa-pf" (if you use partition function) or "-aligner=picxaa-phmm" (if you use pair HMM in order to evaluate the alignment).
          * Divided "-trimqual" parameter (to trim nucleotides according to quality) as "-trimqualleft=(Number)" and "-trimqualright=(Number)".
          * Fixed minor error about colon in FASTA files after processing in PrinSeq.
          * "-design=" and "-br[y|n]" are only mandatory if you want to split your dataset into samples. If not, you MUST to write "-nosplit" parameter in the batch mode.

v 1.0.0 - Original version as a Perl script.

REFERENCES:

	* Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ (1990) Basic local alignment search tool. Journal of Molecular Biology 215: 403-10.
	* Bradley RK, Roberts A, Smoot M, Juvekar S, Do J, Dewey C, Holmes I, Pachter L (2009) Fast Statistical Alignment. PLoS Computational Biology 5: e1000392.
	* Charif D, Lobry JR (2007) SeqinR 1.0-2: a contributed package to the R project for statistical computing devoted to biological sequences retrieval and analysis. IN Bastolla U, Porto M, Roman HE, Vendruscolo M. Structural approaches to sequence evolution: Molecules networks and analysis. Springer-Verlag, New York (USA), 207-32.
	* Clement M, Posada D, Crandall K (2000) TCS: a computer program to estimate gene genealogies. Molecular Ecology 9: 1657-1660 
	* Cock PJA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B, de Hoon MJL (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 25(11): 1422-23.
	* Edgar RC (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic Acids Research 32(5): 1792-7.
	* Edgar RC (2010) Search and clustering orders of magnitude faster than BLAST. Bioinformatics 26(19): 2460-1.
	* Edgar RC, Haas BJ, Clemente JC, Quince C, Knight R (2011) UCHIME improves sensitivity and speed of chimera detection. Bioinformatics 27(16): 2194-200.
	* Finn RD, Clements J, Eddy SR (2011) HMMER web server: interactive sequence similarity searching. Nucleic Acids Research 39: W29-W37.
	* Huse SM, Huber JA, Morrison HG, Sogin ML, Welch DM (2007) Accuracy and quality of massively parallel DNA pyrosequencing. Genome Biology 8: R143.
	* Katoh K, Kuma K, Toh H, Miyata T (2005) MAFFT version 5: improvement in accuracy of multiple sequence alignment. Nucleic Acids Research 33: 511-18.
	* Katoh K, Misawa K, Kuma K, Miyata T (2002) MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform. Nucleic Acids Research 30: 3059-66.
	* Katoh K, Standley DM (2013) MAFFT multiple sequence alignment software Version 7: Improvements in performance and usability. Molecular Biology and Evolution 30(4): 772-80.
	* Katoh K, Toh H (2007) PartTree: an algorithm to build an approximate tree from a large number of unaligned sequences. Bioinformatics 23(3):372-4.
	* Lassmann T, Sonnhammer ELL (2005) Kalign - an accurate and fast multiple sequence alignment algorithm. BMC Bioinformatics 6: 298.
	* Li W, Fu L, Niu B, Wu S, Wooley J (2012) Ultrafast clustering algorithms for metagenomic sequence analysis. Briefings in Bioinformatics 13(6): 656-68.
	* Liu Y, Schmidt B, Maskell DL (2010) MSAProbs: multiple sequence alignment based on pair hidden Markoy models and partition function posterior probabilities. Bioinformatics 26(16):1958-64.
	* Loytynoja A, Goldman N (2005) An algorithm for progressive multiple alignment of sequences with insertions. Proceedings of the National Academy of Sciences of USA 102: 10557-62.
	* Lyras DP, Metzler (2014) ReformAlign: Improved multiple sequence alignments using a profile-based meta-alignment approach". BMC Bioinformatics 15:265. 
	* Martin M (2011) Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal 17.
	* Morgulis A, Gertz EM, Schaffer AA, Agarwala R (2006) A fast and symmetric DUST implementation to mask low-complexity DNA sequences. Journal of Computational Biology 13: 1028-40.
	* Oliphant TE (2007) Python for scientific computing. Computing in Science & Engineering 9: 10-20.
	* Paradis E, Claude J, Strimmer K (2004) APE: analyses of phylogenetics and evolution in R language. Bioinformatics 20: 289-290.
	* Russell DJ, Otu HH, Sayood K (2008) Grammar-based distance in progressive multiple sequence alignment. BMC Bioinformatics 9: 306.
	* Sahraeian SME, Yoon BJ (2010) PicXAA: greedy probabilistic construction of maximum expected accuracy alignment of multiple sequences. Nucleic Acids Research 38(15): 4917-28.
	* Schmieder R, Edwards R (2011) Quality control and preprocessing of metagenomic datasets. Bioinformatics 27: 863-4.
	* Schwartz AS, Pachter L (2007) Multiple alignment by sequence annealing. Bioinformatics 23: e24-9.
	* Sievers F, Dineen D, Wilm A, Higgins DG (2013). Making automated multiple alignments of very large numbers of protein sequences. Bioinformatics 29(8): 989-95.
	* Sievers F, Wilm A, Dineen DG, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H, Remmert M, Soeding J, Thompson JD, Higgins DG (2011). Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. Molecular Systems Biology 7: 539.
	* Templeton AR, Crandall KA, Sing CF (1992) A cladistic analysis of phenotypic associations with haplotypes inferred from restriction endonuclease mapping and DNA sequence data. III. Cladogram estimation. Genetics 132: 619-33.
	* Wheeler TJ, Kececioglu JD (2007) Multiple alignment by aligning alignments. Bioinformatics 23: i559-i568.
	* Wirawan A, Harris RS, Liu Y, Schmidt B, Schröder J (2014) HECTOR: a parallel multistage homopolymer spectrum based error corrector for 454 sequencing data. BMC Bioinformatics 2014, 15:131.
	* Zhang J, Kobert K, Flouri T, Stamatakis A (2014) PEAR: a fast and accurate Illumina Paired-End reAd mergeR. Bioinformatics 30(5): 614-20.
''')

## Filtering HMMER results function
def HMMfiltering_SingleHMM(hmmtable, fastafile, hmmcutoff, paired):
	
	# Naming the new file
	newfastafile = fastafile.replace(".twice.fasta",".hmmfiltered.fasta")
	
	# Starting variables and dictionaries
	evalue = {}
	hmmcutoff = float(hmmcutoff)
	lines_seen = set()
	
	# Saving the list into the memory
	my_hmmtable = open(hmmtable, "rU") if os.path.isfile(hmmtable) else sys.exit("There is no %s in your files\n" % hmmtable)
	infotable = my_hmmtable.read().splitlines()
	for line in infotable:
		if line not in lines_seen:
			if re.match(r"^[^#]", line):
				lines_seen.add(line)
				tbline = re.split(r"\s+", line)
				if tbline[8] == "+":
					evalue[tbline[0]] = float(tbline[4])
	my_hmmtable.close()

	# Saving the HMM-filtered fasta file
	origfasta = open(fastafile, "rU") if os.path.isfile(fastafile) else sys.exit("%s can not be found\n" % fastafile)
	newfasta_handle = open(newfastafile, "w")
	for record in SeqIO.parse(origfasta, "fasta"):
		try:
			if evalue[record.id] < hmmcutoff:
				SeqIO.write(record, newfasta_handle, "fasta")
		except KeyError: 
			pass
	origfasta.close()
	newfasta_handle.close()

## Modifying headers in original fasta file
def ModifyingHeaders(inputfastafile):
	
	if os.path.isfile(inputfastafile):

		# Naming the output fasta file
		outputfasta = inputfastafile.replace(".fasta", ".mod.fasta")

		# Starting changing headers to be functional in the pipeline
		origfasta = open(inputfastafile, "r")
		modfasta = open(outputfasta, "w")
		for record in SeqIO.parse(origfasta, "fasta"):
			record.id = re.sub("[.:]","_",record.id)
			record.description = record.id
			SeqIO.write(record, modfasta, "fasta")
		origfasta.close()
		modfasta.close()

		# Replacing the old file for the new one
		os.rename(outputfasta, inputfastafile)

	else:
		sys.exit("The program cannot open %s\n" % inputfastafile)

# Sorting using natural numbers
def atoi(text):
	return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split('(\d+)', text) ]

# NJ processing results subroutine
def ProcessingNJresults(inputfastafile, samplesfile, derepfile, namefile):
	
	if os.path.isfile(inputfastafile):

		# Starting new dictionaries
		ClusAbund = {}
		ids = {}
		oldsequencesheaders = {}
		oldsequencesheaders2 = {}
		newsequencesheaders = {}
		samples = {}

		# Naming the output fasta file
		outputfasta = inputfastafile.replace(".fasta", ".mod.fasta")

		# Changing headers to the representative sequences variants FASTA file
		i = 0
		idsample = {}
		occurence = {}
		origfasta = open(inputfastafile, "r")
		modfasta = open(outputfasta, "w")
		for record in SeqIO.parse(origfasta, "fasta"):
			i += 1
			record.id = re.search("(\S+);size=\d+;", record.id)
			idsample[i] = record.id.group(1)
			record.id = "Cluster_%s" % i
			record.description = "Cluster_%s" % i
			SeqIO.write(record, modfasta, "fasta")
		origfasta.close()
		modfasta.close()

		# Replacing the old file for the new one
		os.rename(outputfasta, inputfastafile)

		# Saving the DEREP file in memory
		my_derepfile = open(derepfile) if os.path.isfile(derepfile) else sys.exit("%s can not be found\n" % derepfile)
		derepline = my_derepfile.read().splitlines()
		for line in derepline:
			mainderepline = re.split(r"\t", line)
			refName = mainderepline[0]
			Synonims = mainderepline[1]
			oldsequencesheaders[refName] = Synonims
		my_derepfile.close()

		# Saving the NAME file in memory
		my_namefile = open(namefile) if os.path.isfile(namefile) else sys.exit("%s can not be found\n" % namefile)
		nameline = my_namefile.read().splitlines()
		for line in nameline:
			mainnameline = re.split(r"\t", line)
			SRefName = re.sub(';size=\d+;','',mainnameline[0])
			anotherelements = re.sub(';size=\d+;','',mainnameline[1])
			newsequencesheaders[SRefName] = anotherelements
		my_namefile.close()

		# Saving the SAMPLES file in memory
		my_samples = open(samplesfile) if os.path.isfile(samplesfile) else sys.exit("%s can not be found\n" % samplesfile)
		my_samplesfile = my_samples.read().splitlines()
		for sampleline in my_samplesfile:
			if re.match(r"^[^#]", sampleline):
				mainline = re.split(r"\s+", sampleline)
				sample = mainline[0]
				samples[sample] = 1
		my_samples.close()

		for sample in sorted(samples.keys()):
			ids[sample] = {}
			for network in sorted(idsample.keys()):
				ids[sample][network] = 0

		# Adding all synonims in the networks dictionary
		for i in sorted(idsample.keys()):
			speciallist = idsample[i]
			for k1 in sorted(newsequencesheaders.keys()):
				if str(speciallist) == str(k1):
					speciallist = newsequencesheaders[k1].rsplit(",")
			for j1 in range(0, len(speciallist)):
				for k2 in sorted(oldsequencesheaders.keys()):
					if str(speciallist[j1]) == str(k2):
						speciallist[j1] = oldsequencesheaders[k2].rsplit(",")
			idsample[i] = str(speciallist).replace("[","").replace("]","").replace("'","").replace(" ","")

		# Counting all sequences for each cluster and sample and calculating the number of sequences for each cluster
		for i in sorted(idsample.keys()):
			elements = idsample[i].rsplit(",")
			ClusAbund[i] = len(elements)
			for j in range(0, len(elements)):
				oldelement = re.search("(\S+)_\S+", elements[j])
				newelement = oldelement.group(1)
				for sample in samples.keys():
					if newelement == sample:
						ids[sample][i] += 1

		# Printing the absolute matrix
		absresultsfile = inputfastafile.replace(".fasta", ".abs.csv")
		with open(absresultsfile, "w") as out:
			out.write("Sample\\Cluster")
			for net in sorted(idsample.keys()):
				out.write("\t%i" % net)
			out.write("\n")
			for sample in sorted(ids.keys()):
				out.write("%s" % sample)
				for net in sorted(idsample.keys()):
					out.write("\t%s" % ids[sample][net]) if ids[sample][net] != None else out.write("\t0")
				out.write("\n")

		# Printing the summarised output of these analysis
		abslogfile = inputfastafile.replace(".fasta", ".log.txt")
		with open(abslogfile, "w") as out:
			out.write("\tCluster_ID\tN_Seqs\n")
			for i in sorted(ClusAbund.keys()):
				out.write("Cluster_%s\t%s\n" % (i, ClusAbund[i]))

	else:
		sys.exit("The program cannot open %s\n" % inputfastafile)

## Splitting samples and removing barcodes function
def RemovingPrimers(fastafile, oligosfile, splitsamples, designfile, barcodereverse, mismatches):

	# Starting new dictionaries
	barcode = {}
	constructedfwdprimers = {}
	constructedrvsprimers = {}
	rcprimer = {}
	sample = {}

	# Saving a new variable
	mismatches = int(mismatches)
	listmids = []
	seqadapfor = ""
	seqadaprev = ""
	fwd = ""
	rvs = ""

	# Saving the OLIGOS file in memory
	my_oligos = open(oligosfile) if os.path.isfile(oligosfile) else sys.exit("There is no %s in your files\n" % oligosfile)
	oligosline = my_oligos.read().splitlines()
	for line in oligosline:
		if re.match(r"^[^#]", line):
			mainline = re.split(r"\s+", line)
			if re.match(r"^seqadapfor\t\S+", line):
				seqadapfor = mainline[1]
			elif re.match(r"^seqadaprev\t\S+", line):
				seqadaprev = mainline[1]
			elif re.match(r"^forward\t\S+", line):
				fwd = mainline[1]
			elif re.match(r"^reverse\t\S+", line):
				rvs = mainline[1]
			elif re.match(r"^barcode\t\S+\t\S+", line):
				barname = mainline[2]
				barcode[barname] = mainline[1]
	my_oligos.close()

	if splitsamples == 1:

		# Setting defaults values
		if barcodereverse == '-brn':
			bcdrvs = 0
		elif barcodereverse == '-bry':
			bcdrvs = 1

		# Saving the DESIGN file in memory
		my_design = open(designfile) if os.path.isfile(designfile) else sys.exit("There is no %s in your files\n" % designfile)
		designline = my_design.read().splitlines()
		for line in designline:
			if re.match(r"^[^#]", line):
				mainline = re.split(r"\s+", line)
				if bcdrvs == 0:
					midfwd = mainline[0]
					sample[midfwd] = mainline[1]
				else:
					midfwd = mainline[0]
					midrvs = mainline[1]
					if not midfwd in sample.keys():
						sample[midfwd] = {}
					else:
						pass
					sample[midfwd][midrvs] = mainline[2]
		my_design.close()

		# Reconstruction of primers
		for mid in sorted(barcode.keys()):
			listmids.extend([mid])
			if seqadapfor != None or seqadapfor == "":
				constructedfwdprimers[mid] = seqadapfor + barcode[mid] + fwd
			else:
				constructedfwdprimers[mid] = barcode[mid] + fwd
			if seqadaprev != None or seqadaprev == "":
				if bcdrvs == 1:
					constructedrvsprimers[mid] = seqadaprev + barcode[mid] + rvs
				else:
					constructedrvsprimers[mid] = seqadapfor + barcode[mid] + rvs
			else:
				if bcdrvs == 1:
					constructedrvsprimers[mid] = barcode[mid] + rvs
				else:
					constructedrvsprimers[mid] = rvs
			rcprimer[mid] = str(Seq(constructedrvsprimers[mid]).reverse_complement())

		# Processing FASTA file
		number = None
		if bcdrvs == 1:
			for mid in listmids:
				for dim in listmids:
					try: 
						newfastafile = "%s.fasta" % sample[mid][dim]
						origfasta = open(fastafile, "rU") if os.path.isfile(fastafile) else sys.exit("%s can not be found\n" % fastafile)
						tempfile_handle = open(newfastafile, "w")
						number = 0
						for record in SeqIO.parse(origfasta, "fasta"):
							if re.match("^%s" % constructedfwdprimers[mid], str(record.seq)) and re.search("%s$" % rcprimer[dim], str(record.seq)):
								number = number + 1
								record.id = "%s_%i" % (sample[mid][dim], number)
								SeqIO.write(record, tempfile_handle, "fasta")
							elif mismatches != None:
								try:
									condstat = "agrep('%s', '%s', max = list(a=%i))" % (barcode[mid], str(record.seq[:(len(barcode[mid])+mismatches)]), mismatches)
									condstat2 = "agrep('%s', '%s', max = list(i=%i, d=%i, s=%i))" % (str(Seq(barcode[dim]).reverse_complement()), str(record.seq[-(len(barcode[dim])+mismatches):]), mismatches, mismatches+1, mismatches)
									if r(condstat) and r(condstat2):
										number = number + 1
										record.id = "%s_%i" % (sample[mid][dim], number)
										SeqIO.write(record, tempfile_handle, "fasta")
								except IndexError:
									pass
						tempfile_handle.close()
						origfasta.close()
					except KeyError:
						pass
		else:
			for mid in listmids:
				newfastafile = "%s.fasta" % sample[mid]
				origfasta = open(fastafile, "rU") if os.path.isfile(fastafile) else sys.exit("%s can not be found\n" % fastafile)
				tempfile_handle = open(newfastafile, "w")
				number = 0
				for record in SeqIO.parse(origfasta, "fasta"):
					if re.match("^%s" % constructedfwdprimers[mid], str(record.seq)) and re.search("%s$" % rcprimer[mid], str(record.seq)):
						number = number + 1
						record.id = "%s_%i" % (sample[mid], number)
						SeqIO.write(record, tempfile_handle, "fasta")
					elif mismatches != None:
						try:
							condstat = "agrep('%s', '%s', max = list(a=%i))" % (barcode[mid], str(record.seq[:(len(barcode[mid])+mismatches)]), mismatches)
							condstat2 = "agrep('%s', '%s', max = list(i=%i, d=%i, s=%i))" % (str(Seq(barcode[mid]).reverse_complement()), str(record.seq[-(len(barcode[mid])+mismatches):]), mismatches, mismatches+1, mismatches)
							if r(condstat) and r(condstat2):
								number = number + 1
								record.id = "%s_%i" % (sample[mid], number)
								SeqIO.write(record, tempfile_handle, "fasta")
						except IndexError:
							pass
				tempfile_handle.close()
				origfasta.close()

		# Removing complete primer of each sequence
		if bcdrvs == 1:
			for mid in sorted(constructedfwdprimers.keys()):
				for dim in sorted(rcprimer.keys()):
					try:
						tempfastafile = sample[mid][dim] + ".fasta"
						partfastafile = sample[mid][dim] + ".part.fasta"
						newfastafile = sample[mid][dim] + ".trimmed.fasta"
						save1 = open(partfastafile, "w")
						hellout = open("/dev/null", "w")
						subprocess.call([cutadaptprogram, "-g", "^%s" % constructedfwdprimers[mid], tempfastafile], stdout=save1, stderr=hellout)
						save1.close()
						save2 = open(newfastafile, "w")
						subprocess.call([cutadaptprogram, "-b", rcprimer[dim], partfastafile], stdout=save2, stderr=hellout)
						save2.close()
						hellout.close()
##						os.rename(tempfastafile, newfastafile) # ONLY FOR DEBUGGING
						os.remove(tempfastafile)
						os.remove(partfastafile)
					except KeyError:
						pass
		else:
			for mid in sorted(constructedfwdprimers.keys()):
				try:
					tempfastafile = sample[mid] + ".fasta"
					partfastafile = sample[mid] + ".part.fasta"
					newfastafile = sample[mid] + ".trimmed.fasta"
					save1 = open(partfastafile, "w")
					hellout = open("/dev/null", "w")
					subprocess.call([cutadaptprogram, "-g", "^%s" % constructedfwdprimers[mid], tempfastafile], stdout=save1, stderr=hellout)
					save1.close()
					save2 = open(newfastafile, "w")
					subprocess.call([cutadaptprogram, "-b", rcprimer[mid], partfastafile], stdout=save2, stderr=hellout)
					save2.close()
					hellout.close()
##					os.rename(tempfastafile, newfastafile) # ONLY FOR DEBUGGING
					os.remove(tempfastafile)
					os.remove(partfastafile)
				except KeyError:
					pass

		# Print the number of sequences per each sample
		print("Sample\tSequences")
		if bcdrvs == 1:
			for mid in sorted(constructedfwdprimers.keys()):
				for dim in sorted(rcprimer.keys()):
					try:
						newfastafile = sample[mid][dim] + ".trimmed.fasta"
						nrheaders = CountingHeaders(fastafile=newfastafile)
						print("%s\t%s" % (sample[mid][dim], nrheaders))
					except KeyError:
						pass
		else:
			for mid in sorted(constructedfwdprimers.keys()):
				try:
					newfastafile = sample[mid] + ".trimmed.fasta"
					nrheaders = CountingHeaders(fastafile=newfastafile)
					print("%s\t%s" % (sample[mid][dim], nrheaders))
				except KeyError:
					pass

	else:
		# Removing complete primer of each sequence
		partfastafile = fastafile.replace(".fasta", ".part.fasta")
		newfastafile = fastafile.replace(".fasta", ".trimmed.fasta")
		save1 = open(partfastafile, "w")
		hellout = open("/dev/null", "w")
		subprocess.call([cutadaptprogram, "-g", "^%s" % fwd, fastafile], stdout=save1, stderr=hellout)
		save1.close()
		save2 = open(newfastafile, "w")
		subprocess.call([cutadaptprogram, "-b", rvs, partfastafile], stdout=save2, stderr=hellout)
		save2.close()
		hellout.close()
		os.remove(partfastafile)

## Re-naming headers in original fasta file
def RenamingHeaders(inputfastafile, namesample):
	
	if os.path.isfile(inputfastafile):

		# Naming the output fasta file
		tempoutputfasta = inputfastafile.replace(".fasta", ".mod.fasta")
		outputfasta = tempoutputfasta.replace(".mod.fasta", ".fasta")

		# Starting changing headers to be functional in the pipeline
		origfasta = open(inputfastafile, "r")
		modfasta = open(tempoutputfasta, "w")
		numseqs = 1
		for record in SeqIO.parse(origfasta, "fasta"):
			record.id = "%s_%d" % (namesample, numseqs)
			numseqs += 1
			record.description = record.id
			SeqIO.write(record, modfasta, "fasta")
		origfasta.close()
		modfasta.close()

		# Replacing the old file for the new one
		os.rename(tempoutputfasta, outputfasta)

	else:
		sys.exit("The program cannot open %s\n" % inputfastafile)


## Reverse Complementary sequences function
def RevComSeq(inputfastafile, outputfastafile):

	# Making the Reverse Complementary fasta file
	origfasta = open(inputfastafile, "rU") if os.path.isfile(inputfastafile) else sys.exit("%s can not be found\n" % inputfastafile)
	rcfasta_handle = open(outputfastafile, "w")
	for record in SeqIO.parse(origfasta, "fasta"):
		record.id = "RC_%s" % record.id
		record.description = record.id
		record.seq = record.seq.reverse_complement()
		SeqIO.write(record, rcfasta_handle, "fasta")
	origfasta.close()
	rcfasta_handle.close()

## Sortin sequences function
# Implementation of the original python code from Robert Edgar in order to solve the "-stable" bug.
def Sorting(infastafile, alnfastafile):

	def ReadSeqs(FileName):
		Seqs = {}
		Id = ""
		File = open(FileName)
		while 1:
			Line = File.readline()
			if len(Line) == 0:
				return Seqs
			if len(Line) == 0:
				continue
			if Line[0] == ">":
				Id = Line[1:][:-1]
				Seqs[Id] = ""
			else:
				if Id == "":
					sys.exit("FASTA file '%s' does not start with '>'" % FileName)
				Seqs[Id] = Seqs[Id] + Line[:-1]
		File.close()

	def ReadSeqs2(FileName):
		Seqs = []
		Labels = []
		File = open(FileName)
		while 1:
			Line = File.readline()
			if len(Line) == 0:
				return Labels, Seqs
			Line = Line.strip()
			if len(Line) == 0:
				continue
			if Line[0] == ">":
				Id = Line[1:]
				Labels.append(Id)
				Seqs.append("")
			else:
				i = len(Seqs)-1
				Seqs[i] = Seqs[i] + Line
		File.close()
	
	InLabels, InSeqs = ReadSeqs2(infastafile)
	AlnSeqs = ReadSeqs(alnfastafile)
	tempalnfastafile = "%s_mod" % alnfastafile

	modfasta = open(tempalnfastafile, "w")
	for Label in InLabels:
		if Label not in AlnSeqs.keys():
			sys.exit("Not found in alignment: %s" % Label)
		modfasta.write(">%s" % Label)
		modfasta.write(AlnSeqs[Label])
	modfasta.close()
	os.rename(tempalnfastafile, alnfastafile)
	
## Splitting big fastafiles
def Splitter(fastafile, chunk):

	inputfile = open(fastafile, "rU")
	current_length = 0
	current_file = 1
	outputfile = open("group_%d.fasta" % current_file, "w")
	for record in SeqIO.parse(inputfile, "fasta"):
		if current_length == chunk:
			outputfile.close()
			current_length = 0
			current_file = current_file + 1
			outputfile = open("group_%d.fasta" % current_file, "w")
		SeqIO.write([record], outputfile, "fasta")
		current_length = current_length + 1
	inputfile.close()
	outputfile.close()

## Statistical parsimony script
def StatisticalParsimony(fastafile, limitid, consider_gaps, b, nameoutfile, derepfile, namefile, samplesfile):

	## Algorithm from Templeton et al (1992) Genetics

	# Preliminar functions
	def numeratorintegrand(teoq, j, m, b):
		return teoq * 2 * teoq ** (int(j) - 1) * (1 - teoq) ** (2 * int(m) + 1) * (1 - teoq / int(b)) * (2 - teoq * ((int(b) + 1) / int(b) )) ** (int(j) - 1) * (1 - 2 * teoq * (1 - teoq / int(b)))

	def denominatorintegrand(teoq, j, m, b):
		return 2 * teoq ** (int(j) - 1) * (1 - teoq) ** (2 * int(m) + 1) * (1 - teoq / int(b)) * (2 - teoq * ((int(b) + 1) / int(b))) ** (int(j) - 1) * (1 - 2 * teoq * (1 - teoq / int(b)))

	def prod(iterable):
		return reduce(operator.mul, iterable, 1)

	# Charging required libraries from R and allowing interaction between numpy and Rpy2
	ape = importr('ape', robject_translations = {"delta.plot": "deltaplot", "dist.dna": "dist_DNA", "dist.nodes": "distnodes", "node.depth": "nodedepth", "node.depth.edgelength": "nodedepth_edgelength", "node.height": "nodeheight", "node.height.clado": "nodeheight_clado", "prop.part": "proppart"})
	seqinr = importr('seqinr')
	rpy2.robjects.numpy2ri.activate()

	# Starting all imported functions
	as_matrix = robjects.r['as.matrix']

	# Knowing the alignment length
	origfasta = open(fastafile, "rU") if os.path.isfile(fastafile) else sys.exit("%s can not be found\n" % fastafile)
	for record in AlignIO.parse(origfasta, "fasta"):
		lengthDNA = record.get_alignment_length()
		print("\nAll your alignment have %i bp." % lengthDNA)
	origfasta.close()

	# Saving the names and occurence of these sequences
	NamesID = []
	occurence = {}
	origfasta = open(fastafile, "rU") if os.path.isfile(fastafile) else sys.exit("%s can not be found\n" % fastafile)
	for record in SeqIO.parse(origfasta, "fasta"):
		if re.search("^\S+;size=\d+;", record.id):
			record.id = re.search("(\S+);size=(\d+);", record.id)
			newrecord = record.id.group(1)
			occurence[newrecord] = record.id.group(2)
		if NamesID != None or NamesID != []:
			NamesID.extend([newrecord])
		else:
			NamesID = [newrecord]
	origfasta.close()

	# Calculating distance matrix based on real number of differences
	print("\nFirst step: Calculation of distance matrix\n------------------------------------------")
	alignmentfile = ape.as_DNAbin(seqinr.read_alignment(fastafile, format="fasta", forceToLower="FALSE"))
	mutationmatrix = np.matrix(as_matrix(ape.dist_DNA(alignmentfile, model="N", pairwise_deletion=0)))
	if consider_gaps == 1:
		indelsmatrix = np.matrix(as_matrix(ape.dist_DNA(alignmentfile, model="indel", pairwise_deletion=0)))
		distancematrix = mutationmatrix + indelsmatrix
	else:
		distancematrix = mutationmatrix[:]
	print("\nDONE: Calculated distance matrix for every sequences pair.")

	## Calculate probability of parsimony until a given limit (by default, limitid=0.95).
	Pind = []
	Ppars = 0
	CorPpars = 0
                
	if limitid != None or limitid != "":

		# Calculating parsimony probability until a given limit
		print("\nSecond step: Calculating and pooling all sequences according to parsimony probability\n-------------------------------------------------------------------------------------")
		for j in range(1, lengthDNA):
			m = lengthDNA - j # m is the number of invariant sites, and j is the number of changes in DNA
			A = quad(numeratorintegrand, 0, 1, args = (j, m, float(b)))
			B = quad(denominatorintegrand, 0, 1, args = (j, m, float(b)))
			if Pind != []:
				Pind.extend([1 - A[0] / B[0]])
				Pitmanq.extend([A[0] / B[0]])
			else:
				Pind = [1 - A[0] / B[0]]
				Pitmanq = [A[0] / B[0]]
			if prod(Pind[0:j]) < float(limitid):
				break
			else:
				Ppars = prod(Pind[0:j])
				limitstep = j
		print("\nParsimony probability: %.6f\nNumber of differences: %i" % (Ppars, limitstep))

		# Pooling sequences according to the Uncorrected Parsimony Probability
		Networks = []
		UsedSeqs = []
		for i in range(0, len(distancematrix)):
			if NamesID[i] in UsedSeqs:
				continue
			Networks.append([NamesID[i]])
			UsedSeqs.append(NamesID[i])
			for j in range(0, len(distancematrix) - i):
				if NamesID[j] in UsedSeqs:
					continue
				if distancematrix.item((i,j)) <= limitstep:
					Networks[-1].append(NamesID[j])
					UsedSeqs.append(NamesID[j])

		DNetworks = {}
		for i in range(0, len(Networks)):
			DNetworks[i] = Networks[i]

		# Calculating N (number of sequences in the networks)
		Nnetwork = {}
		for i in sorted(DNetworks.keys()):
			elemsum = []
			for j in range(0, len(DNetworks[i])):
				for k in sorted(occurence.keys()):
					if str(DNetworks[i][j]) == str(k):
						if elemsum != []:
							elemsum.append(int(occurence[k]))
						else:
							elemsum = [int(occurence[k])]
			Nnetwork[i] = sum(elemsum)

		# Searching for the most representative sequence in each network
		MRSN = {}
		for i in sorted(DNetworks.keys()):
			repvalue = 0
			for j in range(0, len(DNetworks[i])):
				for k in sorted(occurence.keys()):					
					if str(DNetworks[i][j]) == str(k):
						currentp = float(int(occurence[k])/int(Nnetwork[i]))
						if currentp > repvalue:
							repvalue = currentp
							MRSN[i] = k

		## Printing a new FASTA file with the most representatives sequences for every network
		MRSNfastafile = nameoutfile + ".SPU.fasta"
		origfasta = open(fastafile, "rU") if os.path.isfile(fastafile) else sys.exit("%s can not be found\n" % fastafile)
		newfasta_handle = open(MRSNfastafile, "w")
		for record in SeqIO.parse(origfasta, "fasta"):
			record.id = re.search("(\S+);size=\d+;", record.id)
			oldrecord = record.id.group(1)
			for i in sorted(MRSN.keys()):
				if oldrecord == MRSN[i]:
					record.id = "Network_%s" % (i+1)
					record.description = "Network_%s" % (i+1)
					SeqIO.write(record, newfasta_handle, "fasta")
		newfasta_handle.close()
		origfasta.close()
		print("\nDONE: Pooled all sequences in different networks according to the parsimony probability.")

		# Calculating representative sequences matrix based on real number of differences
		print("\nThird step: Pooling all representative sequences after correcting the parsimony probability\n-------------------------------------------------------------------------------------------")
		RSNamesID = []
		origfasta = open(MRSNfastafile, "rU") if os.path.isfile(MRSNfastafile) else sys.exit("%s can not be found\n" % MRSNfastafile)
		for record in SeqIO.parse(origfasta, "fasta"):
			if RSNamesID != None or RSNamesID != []:
				RSNamesID.extend([record.id])
			else:
				RSNamesID = [record.id]
		origfasta.close()
		RSalignmentfile = ape.as_DNAbin(seqinr.read_alignment(MRSNfastafile, format="fasta", forceToLower="FALSE"))
		RSmutationmatrix = np.matrix(as_matrix(ape.dist_DNA(RSalignmentfile, model="N", pairwise_deletion=0)))
		if consider_gaps == 1:
			RSindelsmatrix = np.matrix(as_matrix(ape.dist_DNA(RSalignmentfile, model="indel", pairwise_deletion=0)))
			RSdistancematrix = RSmutationmatrix + RSindelsmatrix
		else:
			RSdistancematrix = RSmutationmatrix[:]
		print("\nRe-calculated distance matrix for every representative sequences pair.\n")

		# Calculating non-parsimony probability in order to correct the total number of steps.
		Oddvalue = 0
		for y in range(0, limitstep+1):
			for indexes in itertools.combinations(range(0,limitstep), y):
				newList = Pitmanq[:]
				for i in indexes:
					newList[i] = Pind[i]
				Oddvalue = Oddvalue + prod(newList[0:(len(newList)+1)])
			Pnonpars = 1 - Oddvalue
			if Pnonpars < float(limitid):
				break
			else:
				CorPpars = Pnonpars
				biglimit = limitstep + 1
#				biglimit = limitstep + y # NOTE: This should be used in original algorithm according to Templeton et al. (1992). However, this way could underestimate the number of networks.
		print("'Non-Parsimony' probability: %.6f\n" % CorPpars)

		# Pooling sequences according to the Corrected Parsimony Probability
		RSNetworks = []
		RSUsedSeqs = []
		for i in range(0, len(RSdistancematrix)):
			if RSNamesID[i] in RSUsedSeqs:
				continue
			RSNetworks.append([RSNamesID[i]])
			RSUsedSeqs.append(RSNamesID[i])
			for j in range(0, len(RSdistancematrix) - i):
				if RSNamesID[j] in RSUsedSeqs:
					continue
				if RSdistancematrix.item((i,j)) <= biglimit:
					print("WARNING: %s and %s appears to belong to the same network... Correcting." % (RSNamesID[i], RSNamesID[j]))
					RSNetworks[-1].append(RSNamesID[j])
					RSUsedSeqs.append(RSNamesID[j])

		DRSNetworks = {}
		for i in range(0, len(RSNetworks)):
			DRSNetworks[i] = RSNetworks[i]

		# Recalculating N (number of sequences in the networks)
		NRSnetwork = {}
		for i in sorted(DRSNetworks.keys()):
			elemsum = []
			for j in range(0, len(DRSNetworks[i])):
				Checker = re.search("Network_(\d+)", DRSNetworks[i][j])
				nametemp = int(Checker.group(1)) - 1
				for k in sorted(Nnetwork.keys()):
					if int(nametemp) == int(k):
						if elemsum != []:
							elemsum.append(int(Nnetwork[k]))
						else:
							elemsum = [int(Nnetwork[k])]
			NRSnetwork[i] = sum(elemsum)

		# Searching for the most representative sequence in each network
		MRSN2 = {}
		for i in sorted(DRSNetworks.keys()):
			NRSnetworktemp = 0
			for j in range(0, len(DRSNetworks[i])):
				Checker = re.search("Network_(\d+)", DRSNetworks[i][j])
				nametemp = int(Checker.group(1)) - 1
				if Nnetwork[nametemp] > NRSnetworktemp:
					MRSN2[i] = nametemp
					NRSnetworktemp = Nnetwork[nametemp]
				else:
					continue

		# Printing a new FASTA file with the most representatives sequences for every network
		MRSNfastafile2 = nameoutfile + ".fasta"
		origfasta = open(MRSNfastafile, "rU") if os.path.isfile(MRSNfastafile) else sys.exit("%s can not be found\n" % MRSNfastafile)
		newfasta_handle = open(MRSNfastafile2, "w")
		for record in SeqIO.parse(MRSNfastafile, "fasta"):
			record.id = re.search("Network_(\d+)", record.id)
			oldrecord = int(record.id.group(1)) - 1
			for i in sorted(MRSN2.keys()):
				if oldrecord == MRSN2[i]:
					record.id = "Network_%s" % (i+1)
					record.description = "Network_%s" % (i+1)
					SeqIO.write(record, newfasta_handle, "fasta")
		newfasta_handle.close()
		origfasta.close()
		os.remove(MRSNfastafile)
		print("\nDONE: Pooled all representative sequences after correcting parsimony probability.")

		## Obtaining the absolute matrix
		print("\nFourth step: Summarising all results in an absolute frequencies matrix and a log file\n-------------------------------------------------------------------------------------")

		# Starting new dictionaries
		ids = {}
		latestsequencesheaders = {}
		oldsequencesheaders = {}
		newsequencesheaders = {}
		samples = {}

		# Saving the DEREP file in memory
		my_derepfile = open(derepfile) if os.path.isfile(derepfile) else sys.exit("%s can not be found\n" % derepfile)
		derepline = my_derepfile.read().splitlines()
		for line in derepline:
			mainderepline = re.split(r"\t", line)
			refName = mainderepline[0]
			Synonims = mainderepline[1]
			oldsequencesheaders[refName] = Synonims
		my_derepfile.close()

		# Saving the NAME file in memory
		my_namefile = open(namefile) if os.path.isfile(namefile) else sys.exit("%s can not be found\n" % namefile)
		nameline = my_namefile.read().splitlines()
		for line in nameline:
			mainnameline = re.split(r"\t", line)
			SRefName = re.sub(';size=\d+','',mainnameline[0])
			anotherelements = re.sub(';size=\d+','',mainnameline[1])
			newsequencesheaders[SRefName] = anotherelements
		my_namefile.close()

		# Saving the SAMPLES file in memory
		my_samples = open(samplesfile) if os.path.isfile(samplesfile) else sys.exit("%s can not be found\n" % samplesfile)
		my_samplesfile = my_samples.read().splitlines()
		for sampleline in my_samplesfile:
			if re.match(r"^[^#]", sampleline):
				mainline = re.split(r"\s+", sampleline)
				sample = mainline[0]
				samples[sample] = 1
		my_samples.close()

		for sample in sorted(samples.keys()):
			ids[sample] = {}
			for network in sorted(DRSNetworks.keys()):
#			for network in sorted(DNetworks.keys()): # ONLY FOR DEBUGGING
				ids[sample][network] = 0

		# Synonim table
		for key in newsequencesheaders.keys():
			if re.search(r",", newsequencesheaders[key]):
				marray = newsequencesheaders[key].replace(";","").rsplit(",")
			else:
				marray = [newsequencesheaders[key].replace(";","")]
			key = key.replace(";","")
			latestsequencesheaders[key] = []
			for element in marray:
				element = element.replace(element, oldsequencesheaders[element])
				if re.search(r",", element):				
					templist = element.rsplit(",")
				else:
					templist = [element]
				templist = str(templist).replace("[", "").replace("]", "").replace("'","").replace(" ", "")
				if re.search(r",", templist):
					newarray = templist.rsplit(",") 
				else:
					newarray = [templist]
				latestsequencesheaders[key].extend(newarray)

		# Adding all synonims in the networks dictionary
		for i in sorted(DNetworks.keys()):
			speciallist = DNetworks[i]
			for j, i in enumerate(speciallist):
				for k in sorted(latestsequencesheaders.keys()):					
					if str(speciallist[j]) == str(k):
						speciallist[j] = latestsequencesheaders[k]
			DNetworks[i] = speciallist

		for i in sorted(DRSNetworks.keys()):
			speciallist2 = DRSNetworks[i]
			for j, i in enumerate(speciallist2):
				oldelement = re.search("Network_(\S+)", speciallist2[j])
				newelement = int(oldelement.group(1)) - 1
				for k in DNetworks.keys():
					if newelement == k:
						speciallist2[j] = str(DNetworks[newelement]).replace("[", "").replace("]", "").replace("'","").replace(" ", "")
			DRSNetworks[i] = speciallist2
	
		for i in DRSNetworks.keys():
#		for i in DNetworks.keys(): # ONLY FOR DEBUGGING
			if isinstance(i, int):
				for j in range(0, len(DRSNetworks[i])):
#				for j in range(0, len(DNetworks[i])): # ONLY FOR DEBUGGING
					newelement = DRSNetworks[i][j].split(",")
#					newelement = str(DNetworks[i][j]).split(",") # ONLY FOR DEBUGGING
					DRSNetworks[i][j] = newelement
#					DNetworks[i][j] = newelement # ONLY FOR DEBUGGING

		# Counting all sequences for each network and sample
		for i in DRSNetworks.keys():
#		for i in DNetworks.keys(): # ONLY FOR DEBUGGING
			if isinstance(i, int):
				for j in range(0, len(DRSNetworks[i])):
#				for j in range(0, len(DNetworks[i])): # ONLY FOR DEBUGGING
					elements = DRSNetworks[i][j]
#					elements = DNetworks[i][j] # ONLY FOR DEBUGGING
					for k in range(0, len(elements)):
						oldelement = re.search("(\S+)_\S+", elements[k])
						newelement = oldelement.group(1)
						for sample in samples.keys():
							if newelement == sample:
								ids[sample][i] += 1

		# Printing the absolute matrix
		absresultsfile = nameoutfile + ".abs.csv"
		with open(absresultsfile, "w") as out:
			out.write("Sample\\Network")
			for net in DRSNetworks.keys():
#			for net in DNetworks.keys(): # ONLY FOR DEBUGGING
				if isinstance(net, int):
					out.write("\t%i" % (net+1))
			out.write("\n")
			for sample in sorted(ids.keys()):
				out.write("%s" % sample)
				for net in ids[sample].keys():
					if isinstance(net, int):
						if ids[sample][net] != None:
							out.write("\t%s" % ids[sample][net])
						else:
							out.write("\t0")
				out.write("\n")

	# Printing the summarised output of these analysis
	abslogfile = nameoutfile + ".log.txt"
	with open(abslogfile, "w") as out:
		out.write("Network_ID\tN_Seqs\n")
		for i in sorted(NRSnetwork.keys()):
#		for i in sorted(Nnetwork.keys()): # ONLY FOR DEBUGGING
			out.write("Network_%s\t%s\n" % (i+1, NRSnetwork[i]))
#			out.write("Network_%s\t%s\n" % (i+1, Nnetwork[i])) # ONLY FOR DEBUGGING
	print("\nDONE: Created absolute frequencies matrix and log file.")

## Analyzing USEARCH tables
def USEARCHanalysis(usearchtable, newlistfile):

	#Starting new dictionaries
	abcltr = {}
	listcltr = {}
	nmcltr = {}

	# Analyzing the USEARCH clusters file
	my_utable = open(usearchtable) if os.path.isfile(usearchtable) else sys.exit("%s can not be found\n" % usearchtable)
	utableline = my_utable.read().splitlines()
	for line in utableline:
		if re.match(r"^[^#]", line):
			if re.match(r"^C", line):
				mainline = re.split(r"\s+", line)
				try:
					repseq = mainline[8]
					nmcltr[repseq] = mainline[1]
					abcltr[repseq] = mainline[2]
				except IndexError:
					pass
	my_utable.close()

	# Re-analyzing the USEARCH clusters file
	my_utable_2 = open(usearchtable) if os.path.isfile(usearchtable) else sys.exit("%s can not be found\n" % usearchtable)
	utable2line = my_utable_2.read().splitlines()
	for line in utable2line:
		if re.match(r"^[^#]", line):
			if re.match(r"^S\s+", line):
				mainline = re.split(r"\s+", line)
				try:
					if nmcltr[mainline[8]] != None:
						listcltr[mainline[8]] = mainline[8]
				except KeyError:
					continue
			elif re.match(r"^H", line):
				mainline = re.split(r"\s+", line)
				try:
					newseq = mainline[8]
					if nmcltr[mainline[9]]:
						if listcltr[mainline[9]] != None:
							listcltr[mainline[9]] = listcltr[mainline[9]] + "," + newseq
						else:
							listcltr[mainline[9]] = mainline[9]
				except KeyError:
					continue
	my_utable_2.close()

	# Saving all results to the new file
	with open(newlistfile, "w") as out:
		for key in listcltr.keys():
			output = "%s\t%s\n" % (key, listcltr[key])
			out.write(output)

## Searching for program routine
def which(program):
	def is_exe(fpath):
		return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

	fpath, fname = os.path.split(program)
	if fpath:
		if is_exe(program):
			return program
	else:
		for path in os.environ["PATH"].split(os.pathsep):
			path = path.strip('"')
			exe_file = os.path.join(path, program)
			if is_exe(exe_file):
				return exe_file

	return None

### Main program

## Starting computer-related variables
RunningOS = platform.system()
ncpu = multiprocessing.cpu_count()
totalmemory = int(psutil.virtual_memory().total / (1024 ** 2))
freememory = int(psutil.virtual_memory().free / (1024 ** 2))

## Starting initial variables of the script
directory = None
nrparameters = len(sys.argv)
option = None

## Print header of the program
print('''
        #########################################
	#	      QRS v. 1.4.0       	#
	#---------------------------------------#
	#   Author: Enrique Gonzalez-Tortuero   #
	#########################################

QRS v.1.4  Copyright (C) 2014 Enrique Gonzalez-Tortuero
This program comes with ABSOLUTELY NO WARRANTY; for details type 'y' in the interface or '--GNU3' in batch mode.
This is free software, and you are welcome to redistribute it under certain conditions; see above for details.
Running in %s with %d Mb RAM (%d freememory Mb free) and %d processors.
''' % (RunningOS, totalmemory, freememory, ncpu))

## Implementing several programs to the pipeline

# General programs
cutadaptprogram = which("cutadapt")
hectorprogram = which("hector")
hmmerprogram = which("hmmbuild")
hmmerprogram2 = which("nhmmer")
prinseqprogram1 = which("prinseq-lite.pl")
prinseqprogram2 = which("prinseq-graphs-noPCA.pl")
pearprogram = which("pear")
usearchprogram = which("usearch")

# Aligners programs
clustaloprogram = which("clustalo")
fsaprogram = which("fsa")
gramalignprogram = which("GramAlign")
kalignprogram = which("kalign")
mafftprogram = which("mafft")
muscleprogram = which("muscle")
opalprogram = which("opal")
picxaaprogram = which("picxaa")
prankprogram = which("prank")

# Meta-aligner program
reformalprogram = which("ReformAlign")

# Checking the installation of the programs
if (cutadaptprogram == None or cutadaptprogram == "") or (hectorprogram == None  or hectorprogram == "") or (hmmerprogram == None  or hmmerprogram == "" or hmmerprogram2 == None or hmmerprogram2 == "") or (prinseqprogram1 == None or prinseqprogram2 == None or prinseqprogram1 == "" or prinseqprogram2 == "") or (pearprogram == None or pearprogram == "") or (usearchprogram == None or usearchprogram == ""):
	if cutadaptprogram == None or cutadaptprogram == "":
		sys.exit("ERROR: CutAdapt program is not installed in your computer. This program MUST be installed BEFORE executing QRS pipeline.")
	elif hectorprogram == None or hectorprogram == "":
		sys.exit("ERROR: HECTOR program is not installed in your computer. This program MUST be installed BEFORE executing QRS pipeline.")
	elif hmmerprogram == None or hmmerprogram2 == None or hmmerprogram == "" or hmmerprogram2 == "":
		sys.exit("ERROR: HMMER 3.1b software is not installed in your computer. This program MUST be installed BEFORE executing QRS pipeline.")
	elif prinseqprogram1 == None or prinseqprogram2 == None or prinseqprogram1 == "" or prinseqprogram2 == "":
		sys.exit("ERROR: PrinSeq Perl Scripts are not installed in your computer. This program MUST be add to the PATH as EXECUTABLES files BEFORE executing QRS pipeline.")
	elif pearprogram == None or pearprogram == "":
		sys.exit("ERROR: PEAR program is not installed in your computer. This program MUST be installed BEFORE executing QRS pipeline.")
	elif usearchprogram == None or usearchprogram == "":
		sys.exit("ERROR: USEARCH is not installed in your computer. This program MUST be installed BEFORE executing QRS pipeline.")
elif (clustaloprogram == None or clustaloprogram == "") or (fsaprogram == None or fsaprogram == "") or (gramalignprogram == None or gramalignprogram == "") or (kalignprogram == None or kalignprogram == "") or (mafftprogram == None or mafftprogram == "") or (muscleprogram == None or muscleprogram == "") or (opalprogram == None or opalprogram == "") or (picxaaprogram == None or picxaaprogram == "") or (prankprogram == None or prankprogram == "") or (reformalprogram == None or reformalprogram == ""):
	if clustaloprogram == None or clustaloprogram == "":
		print("WARNING: Clustal Omega program is not installed in your computer. It should be convenient to install to be able to use such software in QRS pipeline. Otherwise, please, don't activate the '-aligner=clustalo' option")
	if fsaprogram == None or fsaprogram == "":
		print("WARNING: FSA program is not installed in your computer. It should be convenient to install to be able to use such software in QRS pipeline. Otherwise, please, don't activate the '-aligner=fsa' option")
	if gramalignprogram == None or gramalignprogram == "":
		print("WARNING: GramAlign program is not installed in your computer. It should be convenient to install to be able to use such software in QRS pipeline. Otherwise, please, don't activate the '-aligner=gramalign' option")
	if kalignprogram == None or kalignprogram == "":
		print("WARNING: KAlign program is not installed in your computer. It should be convenient to install to be able to use such software in QRS pipeline. Otherwise, please, don't activate the '-aligner=kalign' option")
	if mafftprogram == None or mafftprogram == "":
		print("WARNING: MAFFT program is not installed in your computer. It should be convenient to install to be able to use such software in QRS pipeline. Otherwise, please, don't activate the '-aligner=mafft-fast', '-aligner=mafft-ginsi', '-aligner=mafft-linsi', '-aligner=mafft-dpparttree', '-aligner=mafft-fastaparttree' or '-aligner=mafft-parttree' option")
	if muscleprogram == None or muscleprogram == "":
		print("WARNING: MUSCLE program is not installed in your computer. It should be convenient to install to be able to use such software in QRS pipeline. Otherwise, please, don't activate the '-aligner=muscle' or '-aligner=muscle-profile' option")
	if opalprogram == None or opalprogram == "":
		print("WARNING: Opal program is not installed in your computer. It should be convenient to install to be able to use such software in QRS pipeline. Otherwise, please, don't activate the '-aligner=opal' option")
	if picxaaprogram == None or picxaaprogram == "":
		print("WARNING: PicXAA program is not installed in your computer. It should be convenient to install to be able to use such software in QRS pipeline. Otherwise, please, don't activate the '-aligner=picxaa-pf' or '-aligner=picxaa-phmm' option")
	if prankprogram == None or prankprogram == "":
		print("WARNING: PRANK program is not installed in your computer. It should be convenient to install to be able to use such software in QRS pipeline. Otherwise, please, don't activate the '-aligner=prank' option")
	if reformalprogram == None or reformalprogram == "":
		print("WARNING: ReformAlign program is not installed in your computer. It should be convenient to install to be able to use such software in QRS pipeline. Otherwise, please, don't activate the '-reformal=yes' option")

### The first parameter is, always, the step you want to run.
if nrparameters > 1:
	parameter = sys.argv[1]
else:
	parameter = None

### Running the pipeline.
try:
	parameter
except NameError:
	# Interactive mode
	nameofparameters = ["qrs.py"]
	directory = input("In which folder do you want to run all analyses? ").rstrip('\n')
	nameofparameters.extend(["-folder=%s" % directory])

	print('''\nWhat step do you want to run?
	a) Aligning and creating HMM for your reference file/s.
	b) Calculating the basic statistics for your NGS dataset.
	c) Processing all your NGS dataset to obtain an alignment.
	d) Processing your alignment to establish representative sequence variants.
	y) GNU license version 3
	z) Help''')
	letteroption = input("Please, enter a letter: ").rstrip('\n')
	if re.match(r"[aA1]", letteroption):
		option = '-RM'
	elif re.match(r"[bB2]", letteroption):
		option = '-AM1'
	elif re.match(r"[cC3]", letteroption):
		option = '-AM2'
	elif re.match(r"[dD4]", letteroption):
		option = '-AM3'
	elif re.match(r"[yY]", letteroption):
		option = '--GNU3'
	elif re.match(r"[zZ]", letteroption):
		option = '--help'
	else:
		sys.exit("ERROR: The requested option is not present.\n")
	nameofparameters.extend([option])
else:
	nameofparameters = ["qrs.py"]
	if (parameter == None):
		# Interactive mode
		directory = input("In which folder do you want to run all analyses? ").rstrip('\n')
		nameofparameters.extend(["-folder=%s" % directory])

		print('''\nWhat step do you want to run?
	a) Aligning and creating HMM for your reference file/s.
	b) Calculating the basic statistics for your NGS dataset.
	c) Processing all your NGS dataset to obtain an alignment.
	d) Processing your alignment to establish representative sequence variants.
	y) GNU license version 3
	z) Help''')
		letteroption = input("Please, enter a letter: ").rstrip('\n')
		if re.match(r"[aA1]", letteroption):
			option = '-RM'
		elif re.match(r"[bB2]", letteroption):
			option = '-AM1'
		elif re.match(r"[cC3]", letteroption):
			option = '-AM2'
		elif re.match(r"[dD4]", letteroption):
			option = '-AM3'
		elif re.match(r"[yY]", letteroption):
			option = '--GNU3'
		elif re.match(r"[zZ]", letteroption):
			option = '--help'
		else:
			sys.exit("ERROR: The requested option is not present.\n")
		nameofparameters.extend([option])
	else:
		# Batch mode
		i = 1
		while i < nrparameters:
			if re.match(r"^-folder\=", sys.argv[i]):
				directory = sys.argv[i].replace(r'-folder=','')
			elif re.match(r"^-[AR]M\d?", sys.argv[i]) or re.search(r"-h", sys.argv[i]) or re.match(r"--GNU3", sys.argv[i]):
				option = parameter
			i += 1

## Displaying help even is there is no options
if option == None or option == '-h' or option == '--help':
	print("This is the same as running ", " ".join(nameofparameters))
	HelpingU()

## Displaying GNUv3 license
if option == '--GNU3':
	print("This is the same as running ", " ".join(nameofparameters))
	GNUv3()

## Going to the specified directory
direcname = os.getcwd()
if directory == "." or directory == None or directory == "":
	os.chdir(direcname)
else:
	direcname = direcname + "/" + directory
	os.chdir(direcname)

## Making all analyses
if option == '-RM':
	# Searching for the second parameter
	if nrparameters > 2:
		parameter2 = sys.argv[2]
	else:
		parameter2 = None
	
	# Starting variables
	aligner = None
	alignoption = None
	alignswitch = None
	nralignments = None
	reformalswitch = None
	
	# Interactive or batch mode
	try:
		parameter2
	except NameError:
		# Interactive mode
		alignoption = input("Do you want to align all your reference sequences ([yes]/no)? ").strip('\n')
		if re.match(r"[Yy](es|ES)?", alignoption) or alignoption == "" or alignoption == None:
			alignswitch = 1
			aligner = input("What aligner do you want to use (clustalo/fsa/gramalign/kalign/mafft-dpparttree/mafft-fast/mafft-fastaparttree/mafft-ginsi/mafft-linsi/mafft-parttree/muscle/muscle-profile/opal/picxaa-pf/picxaa-phmm/[prank])? ").strip('\n')
			if aligner != None or aligner != "":
				nameofparameters.extend(["-aligner=%s" % aligner])
			else:
				nameofparameters.extend(["-aligner=prank"])
			metalign = input("Do you want to use ReformAlign to post-process the alignment ([yes]/no)? ").strip('\n')
			if re.match(r"[Yy](es|ES)?", metalign) or metalign == None:
				reformalswitch = 1
				nameofparameters.extend(["-reformal=yes"])
			elif re.match(r"[Nn][Oo]?", metalign):
				reformalswitch = 0
				nameofparameters.extend(["-reformal=no"])
		else:
			alignswitch = 0
			nameofparameters.extend(["-noalign"])
		refnam = input("What is/are the name/s of the reference fasta file/s?\nIf there are many reference fasta files, please, put the names of all files separated by a dash (-)\n").strip('\n')
		nameofparameters.extend(["-reffiles=%s" % refnam])
		references = refnam.split('-')
		print("This is the same as running ", " ".join(nameofparameters))
	else:
		if (parameter2 == None):
			# Interactive mode
			alignoption = input("Do you want to align all your reference sequences ([yes]/no)? ").strip('\n')
			if re.match(r"[Yy](es|ES)?", alignoption) or alignoption == "" or alignoption == None:
				alignswitch = 1
				aligner = input("What aligner do you want to use (clustalo/fsa/gramalign/kalign/mafft-dpparttree/mafft-fast/mafft-fastaparttree/mafft-ginsi/mafft-linsi/mafft-parttree/muscle/muscle-profile/opal/picxaa-pf/picxaa-phmm/[prank])? ").strip('\n')
				if aligner != None or aligner != "":
					nameofparameters.extend(["-aligner=%s" % aligner])
				else:
					nameofparameters.extend(["-aligner=prank"])
				metalign = input("Do you want to use ReformAlign to fine-tune the alignment ([yes]/no)? ").strip('\n')
				if re.match(r"[Yy](es|ES)?", metalign) or metalign == None:
					reformalswitch = 1
					nameofparameters.extend(["-reformal=yes"])
				elif re.match(r"[Nn][Oo]?", metalign):
					reformalswitch = 0
					nameofparameters.extend(["-reformal=no"])
			else:
				alignswitch = 0
				nameofparameters.extend(["-noalign"])
			refnam = input("What is/are the name/s of the reference fasta file/s?\nIf there are many reference fasta files, please, put the names of all files separated by a dash (-)\n").strip('\n')
			nameofparameters.extend(["-reffiles=%s" % refnam])
			references = refnam.split('-')
			print("This is the same as running ", " ".join(nameofparameters))
		else:
			# Batch mode
			i = 2
			while i < nrparameters:
				if re.match(r"^-noalign", sys.argv[i]):
					alignswitch = 0
				if re.match(r"^-aligner\=", sys.argv[i]):
					aligner = sys.argv[i].replace(r'-aligner=','')
				if re.match(r"^-reformal\=", sys.argv[i]):
					metalign = sys.argv[i].replace(r'-reformal=','')
					if re.match(r"[Yy](es|ES)?", metalign):
						reformalswitch = 1
					elif re.match(r"[Nn](O|o)?", metalign):
						reformalswitch = 0
				if re.match(r"^-reffiles\=", sys.argv[i]):
					references = sys.argv[i].replace(r'-reffiles=','').split('-')
				i += 1

	# Default values if they are not defined.
	if alignswitch == None or alignswitch == "":
		alignswitch = 1
	if aligner == None or aligner == "":
		aligner = "prank"
	if reformalswitch == None or reformalswitch == "":
		reformalswitch = 1
	nralignments = len(references)
	
	print('''
	################################
	# Creating all reference HMMs. #
	################################
	''')

	#Running the step
	for referencefile in references:

		if alignswitch == 1:

			# Naming all output files
			if re.search(r"\.f(na|as|asta)$", referencefile):
				if re.search(r"\.fna$", referencefile):
					rootname = referencefile.replace(r".fna", "")
				elif re.search(r"\.fas$", referencefile):
					rootname = referencefile.replace(r".fas", "")
				elif re.search(r"\.fasta$", referencefile):
					rootname = referencefile.replace(r".fasta", "")
				alignedreferencefile = "%s.aligned" % rootname
				alignmentfile = "%s.best.fas" % alignedreferencefile
				alnreffile = alignmentfile.replace(r".best.fas",".fasta")
				if (reformalswitch == 1):
					metalnreffile = alnreffile.replace(r".fasta",".reformal.fasta")
					hmmreferencefile = metalnreffile.replace(r".aligned.reformal.fasta",".hmm")
				else:
					hmmreferencefile = alnreffile.replace(r".aligned.fasta",".hmm")
			else:
				sys.exit("All your reference files MUST be in FASTA format!")

			# Aligning step
			print("\t1. Aligning reference file: %s\n" % referencefile)
			nrheaders = CountingHeaders(fastafile=referencefile)
			if aligner == "clustalo":
				nthreads = "--threads=%d" % ncpu
				subprocess.call([clustaloprogram, "-i", referencefile, "-o", alnreffile, "-v", nthreads, "--seqtype=DNA", "--outfmt=fa"])
			elif aligner == "fsa":
				save = open(alnreffile, "w")
				if nrheaders >= 1000:
					subprocess.call([fsaprogram, "--fast", "--log", "1", referencefile], stdout=save)
				else:
					subprocess.call([fsaprogram, "--log", "0", referencefile], stdout=save)
				save.close()
			elif aligner == "gramalign":
				subprocess.call([gramalignprogram, "-i", referencefile, "-o", alnreffile, "-f", "2", "-F", "1"])
			elif aligner == "kalign":
				subprocess.call([kalignprogram, "-i", referencefile, "-o", alnreffile, "-f", "fasta"])
			elif re.match(r"^mafft", aligner):
				save = open(alnreffile, "w")
				if aligner == "mafft-fast":
					if nrheaders >= 5000:
						subprocess.call([mafftprogram, "--nuc", "--retree", "1", "--ep", "0.123", "--thread", str(ncpu), referencefile], stdout=save)
					else:
						subprocess.call([mafftprogram, "--nuc", "--retree", "2", "--ep", "0.123", "--thread", str(ncpu), referencefile], stdout=save)
				elif aligner == "mafft-ginsi":
					if nrheaders <= 200:
						subprocess.call([mafftprogram, "--nuc", "--maxiterate", "1000", "--globalpair", "--ep", "0.123", "--thread", str(ncpu), referencefile], stdout=save)
					else:
						sys.exit("G-INS-i algorithm does not deal with more than 200 sequences. Please, choose another alternative for your data analysis.")
				elif aligner == "mafft-linsi":
					if nrheaders <= 200:
						subprocess.call([mafftprogram, "--nuc", "--maxiterate", "1000", "--localpair", "--ep", "0.123", "--thread", str(ncpu), referencefile], stdout=save)
					else:
						sys.exit("L-INS-i algorithm does not deal with more than 200 sequences. Please, choose another alternative for your data analysis.")
				elif aligner == "mafft-dpparttree":
					if nrheaders < 10000:
						print("WARNING: This method is not recommended for aligning less than 10000 sequences.")
						subprocess.call([mafftprogram, "--nuc", "--retree", "2", "--maxiterate", "0", "--nofft", "--dpparttree", "--alga", "--partsize", "1000", "--ep", "0.123", "--thread", str(ncpu), referencefile], stdout=save)
					elif nrheaders > 50000:
						subprocess.call([mafftprogram, "--nuc", "--retree", "1", "--dpparttree", "--alga", "--partsize", "1000", "--ep", "0.123", "--thread", str(ncpu), referencefile], stdout=save)
					else:
						subprocess.call([mafftprogram, "--nuc", "--retree", "1", "--maxiterate", "0", "--nofft", "--dpparttree", "--alga", "--partsize", "1000", "--ep", "0.123", "--thread", str(ncpu), referencefile], stdout=save)
				elif aligner == "mafft-parttree":
					if nrheaders < 10000:
						subprocess.call([mafftprogram, "--nuc", "--retree", "2", "--maxiterate", "0", "--nofft", "--parttree", "--alga", "--partsize", "1000", "--ep", "0.123", "--thread", str(ncpu), referencefile], stdout=save)
					elif nrheaders > 50000:
						subprocess.call([mafftprogram, "--nuc", "--retree", "1", "--parttree", "--alga", "--partsize", "1000", "--ep", "0.123", "--thread", str(ncpu), referencefile], stdout=save)
					else:
						subprocess.call([mafftprogram, "--nuc", "--retree", "1", "--maxiterate", "0", "--nofft", "--parttree", "--alga", "--partsize", "1000", "--ep", "0.123", "--thread", str(ncpu), referencefile], stdout=save)
				elif aligner == "mafft-fastaparttree":
					if nrheaders < 10000:
						print("WARNING: This method is not recommended for aligning less than 10000 sequences.")
						subprocess.call([mafftprogram, "--nuc", "--retree", "2", "--maxiterate", "0", "--nofft", "--fastaparttree", "--alga", "--partsize", "1000", "--ep", "0.123", "--thread", str(ncpu), referencefile], stdout=save)
					elif nrheaders > 50000:
						subprocess.call([mafftprogram, "--nuc", "--retree", "1", "--fastaparttree", "--alga", "--partsize", "1000", "--ep", "0.123", "--thread", str(ncpu), referencefile], stdout=save)
					else:
						subprocess.call([mafftprogram, "--nuc", "--retree", "1", "--maxiterate", "0", "--nofft", "--fastaparttree", "--alga", "--partsize", "1000", "--ep", "0.123", "--thread", str(ncpu), referencefile], stdout=save)
				else:
					sys.exit("If you want to use MAFFT, please select one of the following options in '-aligner=' parameter: mafft-dpparttree, mafft-fast, mafft-fastaparttree, mafft-ginsi, mafft-linsi or mafft-parttree")
				save.close()
			elif re.match(r"^muscle", aligner):
				if aligner == "muscle":
					if nrheaders < 5000:
						subprocess.call([muscleprogram, "-in", referencefile, "-out", alnreffile])
					elif nrheaders > 10000:
						subprocess.call([muscleprogram, "-in", referencefile, "-out", alnreffile, "-maxiters", "1", "-diags"])
					else:
						subprocess.call([muscleprogram, "-in", referencefile, "-out", alnreffile, "-maxiters", "2"])
					Sorting(referencefile, alnreffile)
				elif aligner == "muscle-profile":
					os.mkdir("SEQS_SPLIT")
					shutil.copyfile(referencefile, "SEQS_SPLIT/%s" % referencefile)
					direc2 = direcname + "/SEQS_SPLIT"
					os.chdir(direc2)
					Splitter(referencefile, 1)
					print("Your fastafile was divided in %s fasta files\n" % nrheaders) 

					specialfileslist = ""
					corealn = referencefile.replace(r".fasta", ".CORE.fasta")
					i = 1
					while i < 11:
						namefile = "group_%i.fasta" % i
						if specialfileslist == None or specialfileslist == "":
							specialfileslist = namefile
						else:
							specialfileslist = specialfileslist + ",%s" % namefile
						i += 1
					Cat(joinedfile=corealn, listfile=specialfileslist)
					speciallist = specialfileslist.split(",")
					for specialfilename in speciallist:
						fullpath = direc2 + "/%s" % specialfilename
						os.remove(fullpath)
					subprocess.call([muscleprogram, "-in", corealn, "-out", alnreffile, "-quiet"])
					print("Aligned the 10 most abundant sequences of the dataset. This alignment will be used as a reference or core alignment against the rest of the sequences.\n")
			
					for restfile in sorted(glob.glob('group_*.fasta'), key=natural_keys):
						subprocess.call([muscleprogram, "-profile", "-in1", alnreffile, "-in2", restfile, "-out", alnreffile, "-quiet"])
						print("\rAligned %s against the reference or core alignment." % restfile)

					Sorting(referencefile, alnreffile)
					shutil.copyfile(alnreffile, "../%s" % alnreffile)
					os.chdir("..")
					shutil.rmtree(direc2)
			elif aligner == "opal":
				subprocess.call([opalprogram, "--mem", "%sG" % str(freememory), "--in", referencefile, "--out", alnreffile, "--align_method", "mixed"])
			elif re.match(r"^picxaa", aligner):
				save = open(alnreffile, "w")
				if aligner == "picxaa-pf":
					subprocess.call([picxaaprogram, "-PF", "-nuc", "-v", referencefile], stdout=save)
				elif aligner == "picxaa-phmm":
					subprocess.call([picxaaprogram, "-PHMM", "-nuc", "-v", referencefile], stdout=save)
				else:
					sys.exit("If you want to use PicXAA, please select one of the following options in '-aligner=' parameter: picxaa-fm or picxaa-phmm")
			elif aligner == "prank":
				inputparameter = "-d=%s" % referencefile
				outputparameter = "-o=%s" % alignedreferencefile
				if nrheaders3 < 100:
					subprocess.call([prankprogram, "-DNA", inputparameter, outputparameter])
				else:
					subprocess.call([prankprogram, "-DNA", inputparameter, outputparameter, "-nobppa", "-nomafft"])
				tempfile = alignedreferencefile + ".best.fas"
				os.rename(tempfile, alnreffile)

			# Post-processing alignment and creating HMM files
			if reformalswitch == 1:
				print("\n\t1b. Post-processing the original alignment using ReformAlign in: %s" % alnreffile)
				subprocess.call([reformalprogram, "-i", referencefile, "-o", metalnreffile, "-a", alnreffile, "-v"])

				print("\n\t2. Creating the HMM for your alignment: %s\n" % metalnreffile)
				subprocess.call([hmmerprogram, "--cpu", str(ncpu), hmmreferencefile, metalnreffile])
			else:
				print("\n\t2. Creating the HMM for your alignment: %s\n" % alnreffile)
				subprocess.call([hmmerprogram, "--cpu", str(ncpu), hmmreferencefile, alnreffile])

			# Cleaning all intermediate files
			os.remove(alnreffile)
			if reformalswitch == 1:
				os.remove(metalnreffile)

		else:
			# Naming all output files
			if re.search(r"\.f(na|as|asta)$", referencefile):
				if re.search(r"\.fna$", referencefile):
					hmmreferencefile = referencefile.replace(r".fna", ".hmm")
				elif re.search(r"\.fas$", referencefile):
					hmmreferencefile = referencefile.replace(r".fas", ".hmm")
				elif re.search(r"\.fasta$", referencefile):
					hmmreferencefile = referencefile.replace(r".fasta", ".hmm")
			else:
				sys.exit("All your reference files MUST be in FASTA format!")

			print("\n\t1. Creating the HMM for your alignment: %s\n" % referencefile)
			subprocess.call([hmmerprogram, "--cpu", str(ncpu), hmmreferencefile, referencefile])

	# Printing the final statement
	if alignswitch == 1:
		message = "\nIt is done!\nAll alignments were made in "
		references = "References:\n"
		if aligner == "clustalo":
			message = message + "Clustal Omega (Sievers et al 2011). "
			references = references + "\tSievers F, Wilm A, Dineen DG, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H, Remmert M, Söding J, Thompson JD, Higgins DG (2011). Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. Molecular Systems Biology 7:539\n"
		elif aligner == "fsa":
			message = message + "FSA (Bradley et al 2009). "
			references = references + "\tBradley RK, Roberts A, Smoot M, Juvekar S, Do J, Dewey C, Holmes I, Pachter L (2009) Fast Statistical Alignment. PLoS Computational Biology. 5:e1000392.\n"
		elif aligner == "gramalign":
			message = message + "GramAlign (Russell et al 2008). "
			references = references + "\tRussell DJ, Otu HH, Sayood K (2008) Grammar-based distance in progressive multiple sequence alignment. BMC Bioinformatics 9:306.\n"
		elif aligner == "kalign":
			message = message + "KAlign (Lassmann & Sonnhammer 2005). "
			references = references + "\tLassmann T, Sonnhammer ELL (2005) Kalign - an accurate and fast multiple sequence alignment algorithm. BMC Bioinformatics 6:298\n"
		elif aligner == "mafft-fast":
			message = message + "MAFFT (Katoh & Standley 2013) using FFT-NS-2 algorithm (Katoh et al. 2002). "
			references = references + "\tKatoh K, Standley DM (2013) MAFFT multiple sequence alignment software Version 7: Improvements in performance and usability. Molecular Biology and Evolution 30(4):772-80.\n\tKatoh K, Misawa K, Kuma K, Miyata T (2002) MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform. Nucleic Acids Research 30:3059-66.\n"
		elif aligner == "mafft-ginsi":
			message = message + "MAFFT (Katoh & Standley 2013) using G-INS-i algorithm (Katoh et al. 2005). "
			references = references + "\tKatoh K, Standley DM (2013) MAFFT multiple sequence alignment software Version 7: Improvements in performance and usability. Molecular Biology and Evolution 30(4):772-80.\n\tKatoh K, Kuma K, Toh H, Miyata T (2005) MAFFT version 5: improvement in accuracy of multiple sequence alignment. Nucleic Acids Research 33:511-18.\n"
		elif aligner == "mafft-linsi":
			message = message + "MAFFT (Katoh & Standley 2013) using L-INS-i algorithm (Katoh et al. 2005). "
			references = references + "\tKatoh K, Standley DM (2013) MAFFT multiple sequence alignment software Version 7: Improvements in performance and usability. Molecular Biology and Evolution 30(4):772-80.\n\tKatoh K, Kuma K, Toh H, Miyata T (2005) MAFFT version 5: improvement in accuracy of multiple sequence alignment. Nucleic Acids Research 33:511-18.\n"
		elif aligner == "mafft-dpparttree":
			message = message + "MAFFT (Katoh & Standley 2013) with DP-PartTree algorithm (Katoh & Toh 2007). "
			references = references + "\tKatoh K, Standley DM (2013) MAFFT multiple sequence alignment software Version 7: Improvements in performance and usability. Molecular Biology and Evolution 30(4):772-80.\n\tKatoh K, Toh H (2007) PartTree: an algorithm to build an approximate tree from a large number of unaligned sequences. Bioinformatics 23(3):372-4.\n"
		elif aligner == "mafft-parttree":
			message = message + "MAFFT (Katoh & Standley 2013) with PartTree algorithm (Katoh & Toh 2007). "
			references = references + "\tKatoh K, Standley DM (2013) MAFFT multiple sequence alignment software Version 7: Improvements in performance and usability. Molecular Biology and Evolution 30(4):772-80.\n\tKatoh K, Toh H (2007) PartTree: an algorithm to build an approximate tree from a large number of unaligned sequences. Bioinformatics 23(3):372-4.\n"
		elif aligner == "mafft-fastaparttree":
			message = message + "MAFFT (Katoh & Standley 2013) with FastaPartTree algorithm (Katoh & Toh 2007). "
			references = references + "\tKatoh K, Standley DM (2013) MAFFT multiple sequence alignment software Version 7: Improvements in performance and usability. Molecular Biology and Evolution 30(4):772-80.\n\tKatoh K, Toh H (2007) PartTree: an algorithm to build an approximate tree from a large number of unaligned sequences. Bioinformatics 23(3):372-4.\n"
		elif aligner == "muscle":
			message = message + "MUSCLE (Edgar 2004). "
			references = references + "\tEdgar RC (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic Acids Research 32(5):1792-7.\n"
		elif aligner == "muscle-profile":
			message = message + "MUSCLE (Edgar 2004) with the following strategy: the sequences were splitted in two subsets: i) 'Core': it has the ten most abundant sequences and ii) 'Rest'. After that, 'Core' was aligned using the mentioned program. Later, all sequences for 'Rest' subset are sorted by abundance and every sequence is aligned against 'Core' individually. This strategy is recommended in huge datasets due to the high accuracy of the alignment (Sievers et al. 2013). "
			references = references + "\tEdgar RC (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic Acids Research 32(5):1792-7.\n\tSievers F, Dineen D, Wilm A, Higgins DG (2013). Making automated multiple alignments of very large numbers of protein sequences. Bioinformatics 29(8): 989-95.\n"
		elif aligner == "opal":
			message = message + "Opal (Wheeler & Kececioglu 2007). "
			references = references + "\tWheeler TJ, Kececioglu JD (2007) Multiple alignment by aligning alignments. Bioinformatics 23: i559-i568.\n"
		elif aligner == "picxaa-pf" or aligner == "picxaa-phmm":
			message = message + "PicXAA (Sahraeian & Yoon 2010). "
			references = references + "\tSahraeian SME, Yoon BJ (2010) PicXAA: greedy probabilistic construction of maximum expected accuracy alignment of multiple sequences. Nucleic Acids Research 38(15): 4917-28.\n"
		elif aligner == "prank":
			message = message + "PRANK (Löytynoja & Goldman 2008). "
			references = references + "\tLöytynoja A, Goldman N (2008) A model of evolution and structure for multiple sequence alignment. Philosophical Transactions of the Royal Society B 363(1512):3913-9.\n"

		if reformalswitch == 1:
			message = message + "Then, all alignments were automatically post-processed using ReformAlign (Lyras & Metzler 2014). "
			references = references + "\tLyras DP, Metzler D (2014) ReformAlign: Improved multiple sequence alignments using a profile-based meta-alignment approach. BMC Bioinformatics 15:265\n"
		message = message + "Finally, HMM was made using HMMER (Finn et al 2011)."
		references = references + "\tFinn RD, Clements J, Eddy SR (2011) HMMER web server: interactive sequence similarity searching. Nucleic Acids Research 39:W29-W37.\n"
		sys.exit("%s\n%s" % (message, references))
	else:
		message = "\nIt is done! HMM files from all your alignments were made using HMMER (Finn et al 2011).\n"
		references = "References:\n\tFinn RD, Clements J, Eddy SR (2011) HMMER web server: interactive sequence similarity searching. Nucleic Acids Research 39:W29-W37.\n"
		sys.exit("%s\n%s" % (message, references))
	
elif option == '-AM1':
	# Searching for the second parameter
	if nrparameters > 2:
		parameter2 = sys.argv[2]
	else:
		parameter2 = None

	# Starting variables
	informatfile = None
	NGSfastafile = None
	NGSfastqfile = None
	NGSfastq2file = None
	NGSqualityfile = None
	pairedswitch = None
	phred64switch = None

	# Interactive or batch mode
	try:
		parameter2
	except NameError:
		# Interactive mode
		informatfile = input("\nWhat kind of file(/s) do you have (fasta/fastq)? ").strip('\n')
		nameofparameters.extend(["-informat=%s" % informatfile])
		if informatfile == "fasta":
			NGSfastafile = input("What is its name? ").strip('\n')
			nameofparameters.extend(["-fasta=%s" % NGSfastafile])
			qualexists = input("Do you have a quality file (yes/no)? ").strip('\n')
			if re.match(r"[Yy](ES|es)?", qualexists):
				NGSqualityfile = input("What is its name? ").strip('\n')
				nameofparameters.extend(["-quality=%s" % NGSqualityfile])
		elif informatfile == "fastq":
			NGSfastqfile = input("What is its name? ").strip('\n')
			nameofparameters.extend(["-fastq=%s" % NGSfastqfile])
			paired = input("Are they paired (yes/no)? ").strip('\n')
			if re.match(r"[Yy](ES|es)?", paired):
				pairedswitch = 1
				nameofparameters.extend(["-paired=yes"])
				NGSfastq2file = input("What is its name? ").strip('\n')
				nameofparameters.extend(["-fastq2=%s" % NGSfastq2file])
			elif re.match(r"[Nn][Oo]?", paired):
				pairedswitch = 0
				nameofparameters.extend(["-paired=no"])
			else:
				sys.exit("I don't understand your answer.")
			phred64 = input("Have your quality data got Phred+64 format (yes/no)? ").strip('\n')
			if re.match(r"[Yy](ES|es)?", phred64):
				phred64switch = 1
				nameofparameters.extend(["-phred64=yes"])
			else:
				phred64switch = 0
				nameofparameters.extend(["-phred64=no"])
		else:
			sys.exit("This pipeline only works with FASTA or FASTQ files.")
		print("This is the same as running ", " ".join(nameofparameters))

	else:
		if (parameter2 == None):
			# Interactive mode
			informatfile = input("\nWhat kind of file(/s) do you have (fasta/fastq)? ").strip('\n')
			nameofparameters.extend(["-informat=%s" % informatfile])
			if informatfile == "fasta":
				NGSfastafile = input("What is its name? ").strip('\n')
				nameofparameters.extend(["-fasta=%s" % NGSfastafile])
				qualexists = input("Do you have a quality file (yes/no)? ").strip('\n')
				if re.match(r"[Yy](ES|es)?", qualexists):
					NGSqualityfile = input("What is its name? ").strip('\n')
					nameofparameters.extend(["-quality=%s" % NGSqualityfile])
			elif informatfile == "fastq":
				NGSfastqfile = input("What is its name? ").strip('\n')
				nameofparameters.extend(["-fastq=%s" % NGSfastqfile])
				paired = input("Are they paired (yes/no)? ").strip('\n')
				if re.match(r"[Yy](ES|es)?", paired):
					pairedswitch = 1
					nameofparameters.extend(["-paired=yes"])
					NGSfastq2file = input("What is its name? ").strip('\n')
					nameofparameters.extend(["-fastq2=%s" % NGSfastq2file])
				elif re.match(r"[Nn][Oo]?", paired):
					pairedswitch = 0
					nameofparameters.extend(["-paired=no"])
				else:
					sys.exit("I don't understand your answer.")
				phred64 = input("Have your quality data got Phred+64 format (yes/no)? ").strip('\n')
				if re.match(r"[Yy](ES|es)?", phred64):
					phred64switch = 1
					nameofparameters.extend(["-phred64=yes"])
				else:
					phred64switch = 0
					nameofparameters.extend(["-phred64=no"])
			else:
				sys.exit("This pipeline only works with FASTA or FASTQ files.")
			print("This is the same as running ", " ".join(nameofparameters))
		else:
			# Batch mode
			i = 2
			while i < nrparameters:
				if re.match(r"^-informat\=", sys.argv[i]):
					informatfile = sys.argv[i].replace(r'-aligner=','')
				if re.match(r"^-fasta\=", sys.argv[i]):
					NGSfastafile = sys.argv[i].replace(r'-fasta=','')
				if re.match(r"^-quality\=", sys.argv[i]):
					NGSqualityfile = sys.argv[i].replace(r'-quality=','')
				if re.match(r"^-fastq\=", sys.argv[i]):
					NGSfastqfile = sys.argv[i].replace(r'-fastq=','')
				if re.match(r"^-paired\=", sys.argv[i]):
					paired = sys.argv[i].replace(r'-paired=','')
					if re.match(r"[Yy](ES|es)?", paired):
						pairedswitch = 1
					elif re.match(r"[Nn][Oo]?", paired):
						pairedswitch = 0
				if re.match(r"^-fastq2\=", sys.argv[i]):
					NGSfastq2file = sys.argv[i].replace(r'-fastq2=','')
				if re.match(r"^-phred64\=", sys.argv[i]):
					phred64 = sys.argv[i].replace(r'-phred64=','')
					if re.match(r"[Yy](ES|es)?", paired):
						phred64switch = 1
					elif re.match(r"[Nn][Oo]?", paired):
						phred64switch = 0
				i += 1

	# Checking uncoherencies
	if re.match(r"[^fast(a|q)]", informatfile):
		sys.exit("This pipeline only works with FASTA or FASTQ files.")
	elif informatfile == "fasta" and (NGSfastafile == None or NGSfastafile == ""):
		sys.exit("You MUST to write the name of the FASTA file in the parameter '-fasta='")
	elif informatfile == "fastq" and (NGSfastqfile == None or NGSfastqfile == ""):
		sys.exit("You MUST to write the name of the FASTQ file in the parameter '-fastq='")
	elif pairedswitch == 1 and (NGSfastq2file == None or NGSfastq2file == ""):
		sys.exit("You MUST to write the name of the FASTQ file for paired datasets in the parameter '-fastq2='")

	# Setting default parameters
	if pairedswitch == None or pairedswitch == "":
		pairedswitch = 0
	if phred64switch == None or phred64switch == "":
		phred64switch = 0

	# Preparing the namefiles
	if informatfile == 'fasta':
		if re.search(r"\.fna$", NGSfastafile):
			NGSrootnamefile = NGSfastafile.replace(r'.fna','')
		elif re.search(r"\.fas$", NGSfastafile):
			NGSrootnamefile = NGSfastafile.replace(r'.fas','')
		elif re.search(r"\.fasta$", NGSfastafile):
			NGSrootnamefile = NGSfastafile.replace(r'.fasta','')
		NGSgraphdatafile = NGSrootnamefile + ".gd"
	elif informatfile == 'fastq':
		NGSrootnamefile = NGSfastqfile.replace(r'.fastq','')
		NGSgraphdatafile = NGSrootnamefile + ".gd"
			
	print('''
	########################################################################################
	# Obtaining the basic statistics on the NGS dataset you have and the different graphs. #
	########################################################################################
	''')

	# Running the step
	allparamprogram1 = [prinseqprogram1]
	if informatfile == 'fasta':
		allparamprogram1.extend(["-fasta", NGSfastafile])
		if not NGSqualityfile == None and not NGSqualityfile == "":
			allparamprogram1.extend(["-qual", NGSqualityfile])
	elif informatfile == 'fastq':
		allparamprogram1.extend(["-fastq", NGSfastqfile])
		if pairedswitch == 1:
			allparamprogram1.extend(["-fastq2", NGSfastq2file])
		if phred64switch == 1:
			allparamprogram1.extend(["-phred64"])
	allparamprogram1.extend(["-graph_data", NGSgraphdatafile,"-out_good","null","-out_bad","null","-verbose"])

	subprocess.call(allparamprogram1)
	subprocess.call([prinseqprogram2, "-i", NGSgraphdatafile, "-html_all", "-o", NGSrootnamefile])
	sys.exit("Unable to create PRINSEQ_GRAPHS in your current folder.") if os.path.exists("PRINSEQ_GRAPHS") else os.mkdir("PRINSEQ_GRAPHS")
	specialname = "PRINSEQ_GRAPHS/%s" % NGSrootnamefile
	subprocess.call([prinseqprogram2, "-i", NGSgraphdatafile, "-png_all", "-o", specialname])

	# Printing the final statement
	sys.exit("\nIt is done! You have the basic statistics in the HTML file %s.html. This information was analysed using PrinSeq (Schmieder & Edwards 2011).\nReference:\nSchmieder R, Edwards R (2011) Quality control and preprocessing of metagenomic datasets. Bioinformatics 27:863-4.\n" % NGSrootnamefile)

elif option == '-AM2':
	# Searching for the second parameter
	if nrparameters > 2:
		parameter2 = sys.argv[2]
	else:
		parameter2 = None

	# Starting variables
	aligner = None
	alignswitch = None
	allowedmismatches = None
	allowedns = None
	avelength = None
	barcodedrevprim = None
	cutoff = None
	designfile = None
	filtermethod = None
	filterthreshold = None
	hmmreferencefile = None
	hmmthr = None
	informatfile = None
	maximumgc = None
	maximumlength = None
	minclustersize = None
	minimumgc = None
	minimumlength = None
	minimumqual = None
	namesample = None
	NGSfastafile = None
	NGSfastqfile = None
	NGSfastq2file = None
	NGSqualityfile = None
	nosplit = None
	oligosfile = None
	reformalswitch = None
	trimmingqualleft = None
	trimmingqualright = None
	trimmingtails = None
	pairedswitch = None
	phred64switch = None

	try:
		parameter2
	except NameError:
		# Interactive mode
		print("\nFirst step: knowing all parameters to filtering and trimming sequences according to quality.\nBe sure you see before the basic statistics.\n")
		informatfile = input("What kind of file(/s) do you have (fasta/fastq)? ").strip('\n')
		nameofparameters.extend(["-informat=%s" % informatfile])
		if informatfile == "fasta":
			NGSfastafile = input("What is its name? ").strip('\n')
			nameofparameters.extend(["-fasta=%s" % NGSfastafile])
			qualexists = input("Do you have a quality file (yes/no)? ").strip('\n')
			if re.match(r"[Yy](ES|es)?", qualexists):
				NGSqualityfile = input("What is its name? ").strip('\n')
				nameofparameters.extend(["-quality=%s" % NGSqualityfile])
				qualityswitch = input("Do you want to filter sequences according to a minimum quality (yes/no)? ").strip('\n')
				if re.match(r"[Yy](ES|es)?", qualityswitch):
					minimumqual = input("What is the minimum average quality of the sequences to pass the filters? ").strip('\n')
					nameofparameters.extend(["-minqual=%s" % minimumqual])
				qualityswitch2 = input("Do you want to trim sequences according to a quality threshold (yes/no)? ").strip('\n')
				if re.match(r"[Yy](ES|es)?", qualityswitch2):
					qualityswitch3 = input("Do you want to trim sequences in the forward sense of the sequence (yes/no)? ").strip('\n')
					if re.match(r"[Yy](ES|es)?", qualityswitch3):
						trimmingqualleft = input("How many quality do you consider to remove bad nucleotides in the forward sense of the sequence? ").strip('\n')
						nameofparameters.extend(["-trimqualleft=%s" % trimmingqualleft])
					qualityswitch4 = input("Do you want to trim sequences in the reverse sense of the sequence (yes/no)? ").strip('\n')
					if re.match(r"[Yy](ES|es)?", qualityswitch4):
						trimmingqualright = input("How many quality do you consider to remove bad nucleotides in the reverse sense of the sequence? ").strip('\n')
						nameofparameters.extend(["-trimqualright=%s" % trimmingqualright])
		elif informatfile == "fastq":
			NGSfastqfile = input("What is its name? ").strip('\n')
			nameofparameters.extend(["-fastq=%s" % NGSfastqfile])
			paired = input("Are they paired (yes/no)? ").strip('\n')
			if re.match(r"[Yy](ES|es)?", paired):
				pairedswitch = 1
				nameofparameters.extend(["-paired=yes"])
				NGSfastq2file = input("What is its name? ").strip('\n')
				nameofparameters.extend(["-fastq2=%s" % NGSfastq2file])
			elif re.match(r"[Nn][Oo]?", paired):
				pairedswitch = 0
				nameofparameters.extend(["-paired=no"])
			else:
				sys.exit("I don't understand your answer.")
			phred64 = input("Have your quality data got Phred+64 format (yes/no)? ").strip('\n')
			if re.match(r"[Yy](ES|es)?", phred64):
				phred64switch = 1
				nameofparameters.extend(["-phred64=yes"])
			else:
				phred64switch = 0
				nameofparameters.extend(["-phred64=no"])
			qualityswitch = input("Do you want to filter sequences according to a minimum quality (yes/no)? ").strip('\n')
			if re.match(r"[Yy](ES|es)?", qualityswitch):
				minimumqual = input("What is the minimum average quality of the sequences to pass the filters? ").strip('\n')
				nameofparameters.extend(["-minqual=%s" % minimumqual])
			qualityswitch2 = input("Do you want to trim sequences according to a quality threshold (yes/no)? ").strip('\n')
			if re.match(r"[Yy](ES|es)?", qualityswitch2):
				qualityswitch3 = input("Do you want to trim sequences in the forward sense of the sequence (yes/no)? ").strip('\n')
				if re.match(r"[Yy](ES|es)?", qualityswitch3):
					trimmingqualleft = input("How many quality do you consider to remove bad nucleotides in the forward sense of the sequence? ").strip('\n')
					nameofparameters.extend(["-trimqualleft=%s" % trimmingqualleft])
				qualityswitch4 = input("Do you want to trim sequences in the reverse sense of the sequence (yes/no)? ").strip('\n')
				if re.match(r"[Yy](ES|es)?", qualityswitch4):
					trimmingqualright = input("How many quality do you consider to remove bad nucleotides in the reverse sense of the sequence? ").strip('\n')
					nameofparameters.extend(["-trimqualright=%s" % trimmingqualright])
		else:
			sys.exit("This pipeline only works with FASTA or FASTQ files.")
	
		lengthfilterswitch = input("Do you want to filter sequences according to minimum length (yes/no)? ").strip('\n')
		if re.match(r"[Yy](ES|es)?", lengthfilterswitch):
			minimumlength = input("What is the sequences minimum length (in bp) to pass the filters? ").strip('\n')
			nameofparameters.extend(["-minlen=%s" % minimumlength])
		lengthfilterswitch2 = input("Do you want to filter sequences according to maximum length (yes/no)? ").strip('\n')
		if re.match(r"[Yy](ES|es)?", lengthfilterswitch2):
			maximumlength = input("What is the sequences maximum length (in bp) to pass the filters? ").strip('\n')
			nameofparameters.extend(["-maxlen=%s" % maximumlength])
		gcfilterswitch = input("Do you want to filter sequences according to minimum GC content (yes/no)? ").strip('\n')
		if re.match(r"[Yy](ES|es)?", gcfilterswitch):
			minimumgc = input("What is the sequences minimum GC content (in percentage) to pass the filters? ").strip('\n')
			nameofparameters.extend(["-mingc=%s" % minimumgc])
		gcfilterswitch2 = input("Do you want to filter sequences according to maximum GC content (yes/no)? ").strip('\n')
		if re.match(r"[Yy](ES|es)?", gcfilterswitch2):
			maximumgc = input("What is the sequences maximum GC content (in percentage) to pass the filters? ").strip('\n')
			nameofparameters.extend(["-maxgc=%s" % maximumgc])
		allowedns = int(input("How many uncertain nucleotides (Ns) do you allow in your sequences? ").strip('\n'))
		if allowedns != None or allowedns != "":
			nameofparameters.extend(["-allowns=%s" % allowedns])
		ttswitch = input("Do you want to trim homopolynucleotide tails (eg., poly-As) (yes/no)? ").strip('\n')
		if re.match(r"[Yy](ES|es)?", ttswitch):
			trimmingtails = input("How many homopolynucleotide tails do you allow before trimming sequences in their extremes? ").strip('\n')
			if trimmingtails != None or trimmingtails != "":
				nameofparameters.extend(["-trimtails=%s" % trimmingtails])
		homopolymerswitch = input("Do you want to remove homopolymers (yes/no)? ").strip('\n')
		if re.match(r"[Yy](ES|es)?", homopolymerswitch):
			filtermethod = input("Which 'homopolymers removal' method do you want to use (dust/entropy/hector)? ").strip('\n')
			nameofparameters.extend(["-filmet=%s" % filtermethod])
			if filtermethod == "dust" or filtermethod == "entropy":
				filterthreshold = input("Which threshold do you want to use with %s? " % filtermethod).strip('\n')
				nameofparameters.extend(["-filthr=%s" % filterthreshold])

		print("\nSecond step: knowing the HMM filename to filtering sequences according to a previous reference.\n")
		hmmreferencefile = input("What is the reference HMM file or database name? ").strip('\n')
		nameofparameters.extend(["-hmmfile=%s" % hmmreferencefile])
		hmmthr = input("What is the maximum allowed e-value [1e-10: 10**-10]? If you want to change it, please write as '10**(number)' ").strip('\n')
		if hmmthr != None or hmmthr != "":
			nameofparameters.extend(["-hmmthr=%s" % hmmthr])
	
		print("\nThird step: knowing all parameters and the OLIGOS and DESIGN files to splitting sequences according to the samples and remove primers.\n")
		oligosfile = input("This program needs the information of the primers to remove it in the NGS dataset. It is in a file called OLIGOS with has the following structure:\nseqadapfor\tSEQUENCE\nseqadaprev\tSEQUENCE\nforward\tSEQUENCE\nreverse\tSEQUENCE\nbarcode\tSEQUENCE\tBARID\n'forward' and 'reverse' are mandatory in the file.\nWhat is the OLIGOS file name? ").strip('\n')
		nameofparameters.extend(["-oligos=%s" % oligosfile])
		splitsamplesswitch = input("Do you want to split sequences according to your samples (yes/no)? ").strip('\n')
		if re.match(r"[Yy](ES|es)?", splitsamplesswitch):
			nosplit = 0
			allowedmismatches = input("How many mismatches (in bp) do you allow to detect barcodes in your sequences (0 - INFINITE [1])? ").strip('\n')
			if allowedmismatches != None or allowedmismatches != "":
				nameofparameters.extend(["-allowmis=%s" % allowedmismatches])
			brp = input("Do you have barcoded reverse primers (yes/no)? ").strip('\n')
			if re.match(r"[Yy](ES|es)?", brp):
				barcodedrevprim = "-bry"
			elif re.match(r"[Nn][Oo]?", brp):
				barcodedrevprim = "-brn"
			else:
				sys.exit("I don't understand your answer.")
			nameofparameters.extend([barcodedrevprim])
			designfile = input("Another requered file is DESIGN. It has a table with two or three columns (it depends if you has barcoded reverse primers or not) and as the following structure:\nBARID1\t(BARDID2)\tSAMPLEID\n\nWhat is the DESIGN file name? ").strip('\n')
			nameofparameters.extend(["-design=%s" % designfile])
		else:
			nosplit = 1
			nameofparameters.extend(["-nosplit"])
			namesample = input("What is the name of the sample (this name should be present in the SAMPLES file)? ").strip('\n')
			nameofparameters.extend(["-namesample=%s" % namesample])

		print("\nFourth step: knowing all parameters to denoising your dataset and filtering your dataset.\n")
		cutoff = input("Which cutoff do you want to use to clustering sequences (0.00 - 1.00 [AUTO-DENOISING])? ").strip('\n')
		if cutoff != None or cutoff != "":
			nameofparameters.extend(["-cutoff=%s" % cutoff])
		minclustersize = input("Another QRS's filter is based on cluster's size. We consider that clusters with few sequences are NGS artifacts.\nWhich is the minimum cluster size do you want to allow (by default: 3)? ").strip('\n')
		if minclustersize != None or minclustersize != "":
			nameofparameters.extend(["-minclustersize=%s" % minclustersize])

		print("\nFifth step: Alignment options.\n")
		alignoption = input("Do you want to align all your reference sequences ([yes]/no)? ").strip('\n')
		if re.match(r"[Yy](es|ES)?", alignoption) or alignoption == None:
			alignswitch = 1
			aligner = input("What aligner do you want to use (clustalo/fsa/gramalign/kalign/mafft-dpparttree/mafft-fast/mafft-fastaparttree/mafft-ginsi/mafft-linsi/mafft-parttree/muscle/muscle-profile/opal/picxaa-pf/picxaa-phmm/[prank])? ").strip('\n')
			nameofparameters.extend(["-aligner=%s" % aligner])
			metalign = input("Do you want to use ReformAlign to post-process the alignment ([yes]/no)? ").strip('\n')
			if re.match(r"[Yy](es|ES)?", metalign):
				reformalswitch = 1
				nameofparameters.extend(["-reformal=yes"])
			elif re.match(r"[Nn][Oo]?", metalign):
				reformalswitch = 0
				nameofparameters.extend(["-reformal=no"])
		else:
			alignswitch = 0
			nameofparameters.extend(["-noalign"])

		print("\nThis is the same as running ", " ".join(nameofparameters))
	else:
		if (parameter == None):
			# Interactive mode
			print("\nFirst step: knowing all parameters to filtering and trimming sequences according to quality.\nBe sure you see before the basic statistics.\n")
			informatfile = input("What kind of file(/s) do you have (fasta/fastq)? ").strip('\n')
			nameofparameters.extend(["-informat=%s" % informatfile])
			if informatfile == "fasta":
				NGSfastafile = input("What is its name? ").strip('\n')
				nameofparameters.extend(["-fasta=%s" % NGSfastafile])
				qualexists = input("Do you have a quality file (yes/no)? ").strip('\n')
				if re.match(r"[Yy](ES|es)?", qualexists):
					NGSqualityfile = input("What is its name? ").strip('\n')
					nameofparameters.extend(["-quality=%s" % NGSqualityfile])
					qualityswitch = input("Do you want to filter sequences according to a minimum quality (yes/no)? ").strip('\n')
					if re.match(r"[Yy](ES|es)?", qualityswitch):
						minimumqual = input("What is the minimum average quality of the sequences to pass the filters? ").strip('\n')
						nameofparameters.extend(["-minqual=%s" % minimumqual])
					qualityswitch2 = input("Do you want to trim sequences according to a quality threshold (yes/no)? ").strip('\n')
					if re.match(r"[Yy](ES|es)?", qualityswitch2):
						qualityswitch3 = input("Do you want to trim sequences in the forward sense of the sequence (yes/no)? ").strip('\n')
						if re.match(r"[Yy](ES|es)?", qualityswitch3):
							trimmingqualleft = input("How many quality do you consider to remove bad nucleotides in the forward sense of the sequence? ").strip('\n')
							nameofparameters.extend(["-trimqualleft=%s" % trimmingqualleft])
						qualityswitch4 = input("Do you want to trim sequences in the reverse sense of the sequence (yes/no)? ").strip('\n')
						if re.match(r"[Yy](ES|es)?", qualityswitch4):
							trimmingqualright = input("How many quality do you consider to remove bad nucleotides in the reverse sense of the sequence? ").strip('\n')
							nameofparameters.extend(["-trimqualright=%s" % trimmingqualright])
			elif informatfile == "fastq":
				NGSfastqfile = input("What is its name? ").strip('\n')
				nameofparameters.extend(["-fastq=%s" % NGSfastqfile])
				paired = input("Are they paired (yes/no)? ").strip('\n')
				if re.match(r"[Yy](ES|es)?", paired):
					pairedswitch = 1
					nameofparameters.extend(["-paired=yes"])
					NGSfastq2file = input("What is its name? ").strip('\n')
					nameofparameters.extend(["-fastq2=%s" % NGSfastq2file])
				elif re.match(r"[Nn][Oo]?", paired):
					pairedswitch = 0
					nameofparameters.extend(["-paired=no"])
				else:
					sys.exit("I don't understand your answer.")
				phred64 = input("Have your quality data got Phred+64 format (yes/no)? ").strip('\n')
				if re.match(r"[Yy](ES|es)?", phred64):
					phred64switch = 1
					nameofparameters.extend(["-phred64=yes"])
				else:
					phred64switch = 0
					nameofparameters.extend(["-phred64=no"])
				qualityswitch = input("Do you want to filter sequences according to a minimum quality (yes/no)? ").strip('\n')
				if re.match(r"[Yy](ES|es)?", qualityswitch):
					minimumqual = input("What is the minimum average quality of the sequences to pass the filters? ").strip('\n')
					nameofparameters.extend(["-minqual=%s" % minimumqual])
				qualityswitch2 = input("Do you want to trim sequences according to a quality threshold (yes/no)? ").strip('\n')
				if re.match(r"[Yy](ES|es)?", qualityswitch2):
					qualityswitch3 = input("Do you want to trim sequences in the forward sense of the sequence (yes/no)? ").strip('\n')
					if re.match(r"[Yy](ES|es)?", qualityswitch3):
						trimmingqualleft = input("How many quality do you consider to remove bad nucleotides in the forward sense of the sequence? ").strip('\n')
						nameofparameters.extend(["-trimqualleft=%s" % trimmingqualleft])
					qualityswitch4 = input("Do you want to trim sequences in the reverse sense of the sequence (yes/no)? ").strip('\n')
					if re.match(r"[Yy](ES|es)?", qualityswitch4):
						trimmingqualright = input("How many quality do you consider to remove bad nucleotides in the reverse sense of the sequence? ").strip('\n')
						nameofparameters.extend(["-trimqualright=%s" % trimmingqualright])
			else:
				sys.exit("This pipeline only works with FASTA or FASTQ files.")
	
			lengthfilterswitch = input("Do you want to filter sequences according to minimum length (yes/no)? ").strip('\n')
			if re.match(r"[Yy](ES|es)?", lengthfilterswitch):
				minimumlength = input("What is the sequences minimum length (in bp) to pass the filters? ").strip('\n')
				nameofparameters.extend(["-minlen=%s" % minimumlength])
			lengthfilterswitch2 = input("Do you want to filter sequences according to maximum length (yes/no)? ").strip('\n')
			if re.match(r"[Yy](ES|es)?", lengthfilterswitch2):
				maximumlength = input("What is the sequences maximum length (in bp) to pass the filters? ").strip('\n')
				nameofparameters.extend(["-maxlen=%s" % maximumlength])
			gcfilterswitch = input("Do you want to filter sequences according to minimum GC content (yes/no)? ").strip('\n')
			if re.match(r"[Yy](ES|es)?", gcfilterswitch):
				minimumgc = input("What is the sequences minimum GC content (in percentage) to pass the filters? ").strip('\n')
				nameofparameters.extend(["-mingc=%s" % minimumgc])
			gcfilterswitch2 = input("Do you want to filter sequences according to maximum GC content (yes/no)? ").strip('\n')
			if re.match(r"[Yy](ES|es)?", gcfilterswitch2):
				maximumgc = input("What is the sequences maximum GC content (in percentage) to pass the filters? ").strip('\n')
				nameofparameters.extend(["-maxgc=%s" % maximumgc])
			allowedns = int(input("How many uncertain nucleotides (Ns) do you allow in your sequences? ").strip('\n'))
			if allowedns != None or allowedns != "":
				nameofparameters.extend(["-allowns=%s" % allowedns])
			ttswitch = input("Do you want to trim homopolynucleotide tails (eg., poly-As) (yes/no)? ").strip('\n')
			if re.match(r"[Yy](ES|es)?", ttswitch):
				trimmingtails = input("How many homopolynucleotide tails do you allow before trimming sequences in their extremes? ").strip('\n')
				if trimmingtails != None or trimmingtails != "":
					nameofparameters.extend(["-trimtails=%s" % trimmingtails])
			homopolymerswitch = input("Do you want to remove homopolymers (yes/no)? ").strip('\n')
			if re.match(r"[Yy](ES|es)?", homopolymerswitch):
				filtermethod = input("Which 'homopolymers removal' method do you want to use (dust/entropy/hector)? ").strip('\n')
				nameofparameters.extend(["-filmet=%s" % filtermethod])
				if filtermethod == "dust" or filtermethod == "entropy":
					filterthreshold = input("Which threshold do you want to use with %s? " % filtermethod).strip('\n')
					nameofparameters.extend(["-filthr=%s" % filterthreshold])

			print("\nSecond step: knowing the HMM filename to filtering sequences according to a previous reference.\n")
			hmmreferencefile = input("What is the reference HMM file or database name? ").strip('\n')
			nameofparameters.extend(["-hmmfile=%s" % hmmreferencefile])
			hmmthr = input("What is the maximum allowed e-value [1e-10: 10**-10]? If you want to change it, please write as '10**(number)' ").strip('\n')
			if hmmthr != None or hmmthr != "":
				nameofparameters.extend(["-hmmthr=%s" % hmmthr])
	
			print("\nThird step: knowing all parameters and the OLIGOS and DESIGN files to splitting sequences according to the samples and remove primers.\n")
			oligosfile = input("This program needs the information of the primers to remove it in the NGS dataset. It is in a file called OLIGOS with has the following structure:\nseqadapfor\tSEQUENCE\nseqadaprev\tSEQUENCE\nforward\tSEQUENCE\nreverse\tSEQUENCE\nbarcode\tSEQUENCE\tBARID\n'forward' and 'reverse' are mandatory in the file.\nWhat is the OLIGOS file name? ").strip('\n')
			nameofparameters.extend(["-oligos=%s" % oligosfile])
			splitsamplesswitch = input("Do you want to split sequences according to your samples (yes/no)? ").strip('\n')
			if re.match(r"[Yy](ES|es)?", splitsamplesswitch):
				nosplit = 0
				allowedmismatches = input("How many mismatches (in bp) do you allow to detect barcodes in your sequences (0 - INFINITE [1])? ").strip('\n')
				if allowedmismatches != None or allowedmismatches != "":
					nameofparameters.extend(["-allowmis=%s" % allowedmismatches])
				brp = input("Do you have barcoded reverse primers (yes/no)? ").strip('\n')
				if re.match(r"[Yy](ES|es)?", brp):
					barcodedrevprim = "-bry"
				elif re.match(r"[Nn][Oo]?", brp):
					barcodedrevprim = "-brn"
				else:
					sys.exit("I don't understand your answer.")
				nameofparameters.extend([barcodedrevprim])
				designfile = input("Another requered file is DESIGN. It has a table with two or three columns (it depends if you has barcoded reverse primers or not) and as the following structure:\nBARID1\t(BARDID2)\tSAMPLEID\n\nWhat is the DESIGN file name? ").strip('\n')
				nameofparameters.extend(["-design=%s" % designfile])
			else:
				nosplit = 1
				nameofparameters.extend(["-nosplit"])
				namesample = input("What is the name of the sample (this name should be present in the SAMPLES file)? ").strip('\n')
				nameofparameters.extend(["-namesample=%s" % namesample])


			print("\nFourth step: knowing all parameters to denoising your dataset and filtering your dataset.\n")
			cutoff = input("Which cutoff do you want to use to clustering sequences (0.00 - 1.00 [AUTO-DENOISING])? ").strip('\n')
			if cutoff != None or cutoff != "":
				nameofparameters.extend(["-cutoff=%s" % cutoff])
			minclustersize = input("Another QRS's filter is based on cluster's size. We consider that clusters with few sequences are NGS artifacts.\nWhich is the minimum cluster size do you want to allow (by default: 3)? ").strip('\n')
			if minclustersize != None or minclustersize != "":
				nameofparameters.extend(["-minclustersize=%s" % minclustersize])

			print("\nFifth step: Alignment options.\n")
			alignoption = input("Do you want to align all your reference sequences ([yes]/no)? ").strip('\n')
			if re.match(r"[Yy](es|ES)?", alignoption) or alignoption == None:
				alignswitch = 1
				aligner = input("What aligner do you want to use (clustalo/fsa/gramalign/kalign/mafft-dpparttree/mafft-fast/mafft-fastaparttree/mafft-ginsi/mafft-linsi/mafft-parttree/muscle/muscle-profile/opal/picxaa-pf/picxaa-phmm/[prank])? ").strip('\n')
				nameofparameters.extend(["-aligner=%s" % aligner])
				metalign = input("Do you want to use ReformAlign to post-process the alignment ([yes]/no)? ").strip('\n')
				if re.match(r"[Yy](es|ES)?", metalign):
					reformalswitch = 1
					nameofparameters.extend(["-reformal=yes"])
				elif re.match(r"[Nn][Oo]?", metalign):
					reformalswitch = 0
					nameofparameters.extend(["-reformal=no"])
			else:
				alignswitch = 0
				nameofparameters.extend(["-noalign"])

			print("\nThis is the same as running ", " ".join(nameofparameters))
		else:
			# Batch mode
			i = 2
			while i < nrparameters:
				if re.match(r"^-informat\=", sys.argv[i]):
					informatfile = sys.argv[i].replace(r'-aligner=','')
				if re.match(r"^-fasta\=", sys.argv[i]):
					NGSfastafile = sys.argv[i].replace(r'-fasta=','')
				if re.match(r"^-quality\=", sys.argv[i]):
					NGSqualityfile = sys.argv[i].replace(r'-quality=','')
				if re.match(r"^-fastq\=", sys.argv[i]):
					NGSfastqfile = sys.argv[i].replace(r'-fastq=','')
				if re.match(r"^-paired\=", sys.argv[i]):
					paired = sys.argv[i].replace(r'-paired=','')
					if re.match(r"[Yy](ES|es)?", paired):
						pairedswitch = 1
					elif re.match(r"[Nn][Oo]?", paired):
						pairedswitch = 0
				if re.match(r"^-fastq2\=", sys.argv[i]):
					NGSfastq2file = sys.argv[i].replace(r'-fastq2=','')
				if re.match(r"^-phred64\=", sys.argv[i]):
					phred64 = sys.argv[i].replace(r'-phred64=','')
					if re.match(r"[Yy](ES|es)?", paired):
						phred64switch = 1
					elif re.match(r"[Nn][Oo]?", paired):
						phred64switch = 0
				if re.match(r"^-minlen\=", sys.argv[i]):
					minimumlength = sys.argv[i].replace(r'-minlen=','')
				if re.match(r"^-maxlen\=", sys.argv[i]):
					maximumlength = sys.argv[i].replace(r'-maxlen=','')
				if re.match(r"^-mingc\=", sys.argv[i]):
					minimumgc = sys.argv[i].replace(r'-mingc=','')
				if re.match(r"^-maxgc\=", sys.argv[i]):
					maximumgc = sys.argv[i].replace(r'-maxgc=','')
				if re.match(r"^-minqual\=", sys.argv[i]):
					minimumqual = sys.argv[i].replace(r'-minqual=','')
				if re.match(r"^-trimqualleft\=", sys.argv[i]):
					trimmingqualleft = sys.argv[i].replace(r'-trimqualleft=','')
				if re.match(r"^-trimqualright\=", sys.argv[i]):
					trimmingqualright = sys.argv[i].replace(r'-trimqualright=','')
				if re.match(r"^-trimtails\=", sys.argv[i]):
					trimmingtails = sys.argv[i].replace(r'-trimtails=','')
				if re.match(r"^-allowns\=", sys.argv[i]):
					allowedns = int(sys.argv[i].replace(r'-allowns=',''))
				if re.match(r"^-filmet\=", sys.argv[i]):
					filtermethod = sys.argv[i].replace(r'-filmet=','')
				if re.match(r"^-filthr\=", sys.argv[i]):
					filterthreshold = sys.argv[i].replace(r'-filthr=','')
				if re.match(r"^-hmmfile\=", sys.argv[i]):
					hmmreferencefile = sys.argv[i].replace(r'-hmmfile=','')
				if re.match(r"^-hmmthr\=", sys.argv[i]):
					hmmthr = sys.argv[i].replace(r'-hmmthr=','')
				if sys.argv[i] == "-nosplit":
					nosplit = 1
				if re.match(r"-namesample\=", sys.argv[i]):
					namesample = sys.argv[i].replace(r'-namesample=','')
				if re.match(r"^-oligos\=", sys.argv[i]):
					oligosfile = sys.argv[i].replace(r'-oligos=','')
				if re.match(r"^-design\=", sys.argv[i]):
					designfile = sys.argv[i].replace(r'-design=','')
				if re.match(r"^-br[yn]$", sys.argv[i]):
					barcodedrevprim = sys.argv[i]
				if re.match(r"^-allowmis\=", sys.argv[i]):
					allowedmismatches = int(sys.argv[i].replace(r'-allowmis=',''))
				if re.match(r"^-cutoff\=", sys.argv[i]):
					cutoff = float(sys.argv[i].replace(r'-cutoff=',''))
				if re.match(r"^-minclustersize\=", sys.argv[i]):
					minclustersize = sys.argv[i].replace(r'-minclustersize=','')
				if re.match(r"^-noalign", sys.argv[i]):
					alignswitch = 0
				if re.match(r"^-aligner\=", sys.argv[i]):
					aligner = sys.argv[i].replace(r'-aligner=','')
				if re.match(r"^-reformal\=", sys.argv[i]):
					metalign = sys.argv[i].replace(r'-reformal=','')
					if re.match(r"[Yy](es|ES)?", metalign):
						reformalswitch = 1
					elif re.match(r"[Nn](O|o)?", metalign):
						reformalswitch = 0
					else:
						sys.exit("I don't recognize the ReformAlign parameter")

	# Setting default parameters
	if alignswitch == None or alignswitch == "":
		alignswitch = 1
		if aligner == None or aligner == "":
			aligner = "prank"
		if reformalswitch == None or reformalswitch == "":
			reformalswitch = 1
	if allowedmismatches == None or allowedmismatches == "":
		allowedmismatches = 1
	if hmmthr == None or hmmthr == "":
		hmmthr = 10 ** -10
	if minclustersize == None or minclustersize == "":
		minclustersize = 3
	if nosplit == None or nosplit == "":
		nosplit = 0
	if pairedswitch == None or pairedswitch == "":
		pairedswitch = 0
	if phred64switch == None or phred64switch == "":
		phred64switch = 0


	# Checking uncoherencies
	if re.match(r"[^fast(a|q)]", informatfile):
		sys.exit("This pipeline only works with FASTA or FASTQ files.")
	elif informatfile == "fasta" and (NGSfastafile == None or NGSfastafile == ""):
		sys.exit("You MUST to write the name of the FASTA file in the parameter '-fasta='")
	elif informatfile == "fastq" and (NGSfastqfile == None or NGSfastqfile == ""):
		sys.exit("You MUST to write the name of the FASTQ file in the parameter '-fastq='")
	elif pairedswitch == 1 and (NGSfastq2file == None or NGSfastq2file == ""):
		sys.exit("You MUST to write the name of the FASTQ file for paired datasets in the parameter '-fastq2='")
	elif hmmreferencefile == None or hmmreferencefile == "":
		sys.exit("You MUST to write the name of the HMM REFERENCE file in the parameter '-hmmfile='")
	elif oligosfile == None or oligosfile == "":
		sys.exit("You MUST to write the name of the OLIGOS file in the parameter '-oligos='")
	elif nosplit == 0:
		if designfile == None or designfile == "":
			sys.exit("You MUST to write the name of the DESIGN file in the parameter '-design='")
		elif barcodedrevprim == None or barcodedrevprim == "":
			sys.exit("The presence or absence of barcoded reverse primers was not specified. Please, use the parameters '-bry' or '-brn' if you have (or not) barcoded reverse primers")
	elif nosplit == 1:
		if namesample == None or namesample == "":
			sys.exit("You MUST to define the name of the sample in the parameter '-namesample='if your dataset is separated by samples. This name should be the same as used in the SAMPLES file.")
	
	# Preparing the namefiles
	if informatfile == 'fasta':
		if re.search(r"\.fna$", NGSfastafile):
			NGSrootnamefile = NGSfastafile.replace(r'.fna','')
		elif re.search(r"\.fas$", NGSfastafile):
			NGSrootnamefile = NGSfastafile.replace(r'.fas','')
		elif re.search(r"\.fasta$", NGSfastafile):
			NGSrootnamefile = NGSfastafile.replace(r'.fasta','')
	elif informatfile == 'fastq':
		NGSrootnamefile = NGSfastqfile.replace(r'.fastq','')
	NGSbadsequences = "fakeseqs"
	NGSgoodsequences = "%s.good" % NGSrootnamefile
	NGSlogfile = "%s.txt" % NGSrootnamefile
	NGSrevcompfastafile = "%s.rc.fasta" % NGSgoodsequences
	NGStwicefastafile = "%s.twice.fasta" % NGSgoodsequences
	
	if pairedswitch == 0:
		NGSgoodfastafile = "%s.fasta" % NGSgoodsequences
		NGSgrcfileslist = "%s,%s" % (NGSgoodfastafile, NGSrevcompfastafile)
	else:
		NGSgoodfastqfile1 = "%s_1.fastq" % NGSgoodsequences
		NGSgoodfastqfile2 = "%s_2.fastq" % NGSgoodsequences
		NGSmergedfastqfile = "%s.assembled.fastq" % NGSgoodsequences
		NGSmergedfastafile = "%s.assembled.fasta" % NGSgoodsequences
		NGSgrcfileslist = "%s,%s" % (NGSmergedfastafile, NGSrevcompfastafile)

	HMMresults = NGStwicefastafile.replace(".fasta", ".table")
	NGShmmfilteredfastafile = "%s.hmmfiltered.fasta" % NGSgoodsequences
	tempalltrimmedderepfile = NGShmmfilteredfastafile.replace(r".fasta", ".temptrimmed.fasta")
	alltrimmedfile = NGShmmfilteredfastafile.replace(r".fasta", ".trimmed.fasta")
	
	alltrimmedderepfile = alltrimmedfile.replace(r".fasta", ".derep.fasta")
	alltrimmeddereptable = ("%s.uctbl" % alltrimmedderepfile).replace(r".fasta","")
	alltrimmeddereplist = "derep.list"

	alltrimmednonchimerafile = alltrimmedderepfile.replace(r".fasta", ".nonchimera.fasta")
	alltrimmeduniquetable = ("%s.uctbl" % alltrimmednonchimerafile).replace(r".fasta", "")
	alltrimmedclustlist = "names.list"

	alltrimmedalignment = alltrimmednonchimerafile.replace(r".fasta", ".aligned.fasta")
	metNGSfile = alltrimmedalignment.replace(r".fasta", ".reformal.fasta")

	# Running all steps
	print('''
	###############################################
	# Filtering and aligning all your NGS dataset #
	###############################################
	''')

	# Filtering and trimming sequences
	print("\t1. Filtering and trimming sequences from the NGS data set according to your parameters\n")
	allparamprogram1 = [prinseqprogram1]
	if informatfile == 'fasta':
		allparamprogram1.extend(["-fasta", NGSfastafile])
		if not NGSqualityfile == None and not NGSqualityfile == "":
			allparamprogram1.extend(["-qual", NGSqualityfile])
	elif informatfile == 'fastq':
		allparamprogram1.extend(["-fastq", NGSfastqfile])
		if pairedswitch == 1:
			allparamprogram1.extend(["-fastq2",NGSfastq2file])
		if phred64switch == 1:
			allparamprogram1.extend(["-phred64"])
	if not minimumqual == None and not minimumqual == "":
		allparamprogram1.extend(["-min_qual_mean", minimumqual])
	if not trimmingqualleft == None and not trimmingqualleft == "":
		allparamprogram1.extend(["-trim_qual_left", trimmingqualleft])
	if not trimmingqualright == None and not trimmingqualright == "":
		allparamprogram1.extend(["-trim_qual_right", trimmingqualright])
	if not trimmingtails == None and not trimmingtails == "":
		allparamprogram1.extend(["-trim_tail_left", trimmingtails, "-trim_tail_right", trimmingtails])
	if not allowedns == None and not allowedns == "":
		allparamprogram1.extend(["--ns_max_n", str(allowedns)])
	if not minimumlength == None and not minimumlength == "":
		allparamprogram1.extend(["-min_len", minimumlength])
	if not maximumlength == None and not maximumlength == "":
		allparamprogram1.extend(["-max_len", maximumlength])
	if not minimumgc == None and not minimumgc == "":
		allparamprogram1.extend(["-min_gc", minimumgc])
	if not maximumgc == None and not maximumgc == "":
		allparamprogram1.extend(["-max_gc", maximumgc])
	if (filtermethod == "dust" or filtermethod == "entropy") and (not filterthreshold == None and not filterthreshold == ""):
		allparamprogram1.extend(["-lc_method", filtermethod, "-lc_threshold", filterthreshold])
	if not NGSgoodsequences == None and not NGSgoodsequences == "":
		if pairedswitch == 0:
			allparamprogram1.extend(["-out_format", "1", "-out_good", NGSgoodsequences, "-out_bad", NGSbadsequences, "-log", NGSlogfile, "-seq_case", "upper", "-verbose"])
		else:
			allparamprogram1.extend(["-out_format", "3", "-out_good", NGSgoodsequences, "-out_bad", NGSbadsequences, "-log", NGSlogfile, "-seq_case", "upper", "-verbose"])

	subprocess.call(allparamprogram1)

	# Merging paired end reads using PEAR and/or creating the Reverse Complementary sequences and making a whole dataset to analyse.
	if pairedswitch == 1:
		print("\n\t1b. Merging all pre-filtered NGS sequences from both FASTQ files.\n")
		allparamprogram2 = [pearprogram, "-f", NGSgoodfastqfile1, "-r", NGSgoodfastqfile2, "-o", NGSgoodsequences]
		if phred64switch == 1:
			allparamprogram2.extend(["-b", "64"])
		subprocess.call(allparamprogram2)
		count = SeqIO.convert(NGSmergedfastqfile, "fastq", NGSmergedfastafile, "fasta")
		ModifyingHeaders(inputfastafile=NGSmergedfastafile)
		if filtermethod == "hector":
			print("\n")
			tempfile = NGSmergedfastafile.replace(".fasta", ".temp.fasta")
			subprocess.call([hectorprogram, "-p", str(ncpu), "-o", tempfile, NGSmergedfastafile])
			os.rename(tempfile, NGSmergedfastafile)
		RevComSeq(inputfastafile=NGSmergedfastafile, outputfastafile=NGSrevcompfastafile)
		Cat(joinedfile=NGStwicefastafile, listfile=NGSgrcfileslist)
	else:
		ModifyingHeaders(inputfastafile=NGSgoodfastafile)
		if filtermethod == "hector":
			print("\n")
			tempfile = NGSgoodfastafile.replace(".fasta", ".temp.fasta")
			subprocess.call([hectorprogram, "-p", str(ncpu), "-o", tempfile, NGSgoodfastafile])
			os.rename(tempfile, NGSgoodfastafile)
		RevComSeq(inputfastafile=NGSgoodfastafile, outputfastafile=NGSrevcompfastafile)
		Cat(joinedfile=NGStwicefastafile, listfile=NGSgrcfileslist)

	# Filtering sequences according to the references
	print("\n\t2. Filtering your NGS sequences according to your HMM reference file\n")
	hellout = open("/dev/null", "w")
	subprocess.call([hmmerprogram2, "--toponly", "--dfamtblout", HMMresults, "--cpu", str(ncpu), hmmreferencefile, NGStwicefastafile], stdout=hellout)
	HMMfiltering_SingleHMM(hmmtable=HMMresults, fastafile=NGStwicefastafile, hmmcutoff=hmmthr, paired=pairedswitch)
	hellout.close()
	nrheaders = CountingHeaders(fastafile=NGShmmfilteredfastafile)
	print("In this step, you recovered %d sequences.\n" % nrheaders)

	# Splitting samples and removing barcodes.
	if nosplit == 0:
		print("\t3. Putting sequences into samples and removing primers\n")
		RemovingPrimers(fastafile=NGShmmfilteredfastafile, oligosfile=oligosfile, splitsamples=1, designfile=designfile, barcodereverse=barcodedrevprim, mismatches=str(int(allowedmismatches)))
		fileslist = ""
		for trimmedfile in glob.glob('*.trimmed.fasta'):
			if fileslist == None or fileslist == "":
				fileslist = trimmedfile
			else:
				fileslist = fileslist + ",%s" % trimmedfile
		Cat(joinedfile=alltrimmedfile, listfile=fileslist)
		nrheaders2 = CountingHeaders(fastafile=alltrimmedfile)
		print("In this step, you recovered %d sequences.\n" % nrheaders2)
	else:
		print("\t3. Removing primers using CutAdapt\n")
		RemovingPrimers(fastafile=NGShmmfilteredfastafile, oligosfile=oligosfile, splitsamples=0, designfile=None, barcodereverse=None, mismatches=str(int(allowedmismatches)))
		RenamingHeaders(NGShmmfilteredfastafile, namesample)
		os.rename(NGShmmfilteredfastafile, alltrimmedfile)
		nrheaders2 = CountingHeaders(fastafile=alltrimmedfile)
		print("In this step, you recovered %d sequences.\n" % nrheaders2)

	# Clustering and denoising sequences
	print("\t4. Dereplicating, denoising and removing chimeric sequences from all your sequences\n")
	if cutoff == None or cutoff == "":
		avelength = AveLength(sequencefile=alltrimmedfile)
		cutoff = float((avelength - 1) / avelength - 0.001)
	perid = float(cutoff) * 100
	cutoffperotus = float(100 - perid)
	subprocess.call([usearchprogram, "-derep_fulllength", alltrimmedfile, "-output", tempalltrimmedderepfile, "-uc", alltrimmeddereptable, "-sizeout"])
	subprocess.call([usearchprogram, "-sortbysize", tempalltrimmedderepfile, "-output", alltrimmedderepfile, "-minsize", str(minclustersize), "-quiet"])
	os.remove(tempalltrimmedderepfile)
	USEARCHanalysis(usearchtable=alltrimmeddereptable, newlistfile=alltrimmeddereplist)

	print("\n")
	subprocess.call([usearchprogram, "-cluster_otus", alltrimmedderepfile, "-otus", alltrimmednonchimerafile, "-uc", alltrimmeduniquetable, "-otu_radius_pct", "%.8f" % float(cutoffperotus), "-usersort"])
	USEARCHanalysis(usearchtable=alltrimmeduniquetable, newlistfile=alltrimmedclustlist)

	#Aligning sequences
	if alignswitch == 1:
		print("\n\t5. Aligning file: %s" % alltrimmednonchimerafile)
		nrheaders3 = CountingHeaders(fastafile=alltrimmednonchimerafile)
		if aligner == "clustalo":
			nthreads = "--threads=%d" % ncpu
			subprocess.call([clustaloprogram, "-i", alltrimmednonchimerafile, "-o", alltrimmedalignment, "-v", nthreads, "--seqtype=DNA", "--outfmt=fa"])
		elif aligner == "fsa":
			save = open(alltrimmedalignment, "w")
			if nrheaders3 >= 1000:
				subprocess.call([fsaprogram, "--fast", "--log", "0", alltrimmednonchimerafile], stdout=save)
			else:
				subprocess.call([fsaprogram, "--log", "1", alltrimmednonchimerafile], stdout=save)
			save.close()
		elif aligner == "gramalign":
			subprocess.call([gramalignprogram, "-i", alltrimmednonchimerafile, "-o", alltrimmedalignment, "-f", "2", "-F", "1"])
		elif aligner == "kalign":
			subprocess.call([kalignprogram, "-i", alltrimmednonchimerafile, "-o", alltrimmedalignment, "-f", "fasta"])
		elif re.match(r"^mafft", aligner):
			save = open(alltrimmedalignment, "w")
			if aligner == "mafft-fast":
				if nrheaders3 >= 5000:
					subprocess.call([mafftprogram, "--nuc", "--retree", "1", "--ep", "0.123", "--thread", str(ncpu), alltrimmednonchimerafile], stdout=save)
				else:
					subprocess.call([mafftprogram, "--nuc", "--retree", "2", "--ep", "0.123", "--thread", str(ncpu), alltrimmednonchimerafile], stdout=save)
			elif aligner == "mafft-ginsi":
				if nrheaders3 <= 200:
					subprocess.call([mafftprogram, "--nuc", "--maxiterate", "1000", "--globalpair", "--ep", "0.123", "--thread", str(ncpu), alltrimmednonchimerafile], stdout=save)
				else:
					sys.exit("G-INS-i algorithm does not deal with more than 200 sequences. Please, choose another alternative for your data analysis.")
			elif aligner == "mafft-linsi":
				if nrheaders3 <= 200:
					subprocess.call([mafftprogram, "--nuc", "--maxiterate", "1000", "--localpair", "--ep", "0.123", "--thread", str(ncpu), alltrimmednonchimerafile], stdout=save)
				else:
					sys.exit("L-INS-i algorithm does not deal with more than 200 sequences. Please, choose another alternative for your data analysis.")
			elif aligner == "mafft-dpparttree":
				if nrheaders3 < 10000:
					print("WARNING: This method is not recommended for aligning less than 10000 sequences.")
					subprocess.call([mafftprogram, "--nuc", "--retree", "2", "--maxiterate", "0", "--nofft", "--dpparttree", "--alga", "--partsize", "1000", "--ep", "0.123", "--thread", str(ncpu), alltrimmednonchimerafile], stdout=save)
				elif nrheaders3 > 50000:
					subprocess.call([mafftprogram, "--nuc", "--retree", "1", "--dpparttree", "--alga", "--partsize", "1000", "--ep", "0.123", "--thread", str(ncpu), alltrimmednonchimerafile], stdout=save)
				else:
					subprocess.call([mafftprogram, "--nuc", "--retree", "1", "--maxiterate", "0", "--nofft", "--dpparttree", "--alga", "--partsize", "1000", "--ep", "0.123", "--thread", str(ncpu), alltrimmednonchimerafile], stdout=save)
			elif aligner == "mafft-parttree":
				if nrheaders3 < 10000:
					print("WARNING: This method is not recommended for aligning less than 10000 sequences.")
					subprocess.call([mafftprogram, "--nuc", "--retree", "2", "--maxiterate", "0", "--nofft", "--parttree", "--alga", "--partsize", "1000", "--ep", "0.123", "--thread", str(ncpu), alltrimmednonchimerafile], stdout=save)
				elif nrheaders3 > 50000:
					subprocess.call([mafftprogram, "--nuc", "--retree", "1", "--parttree", "--alga", "--partsize", "1000", "--ep", "0.123", "--thread", str(ncpu), alltrimmednonchimerafile], stdout=save)
				else:
					subprocess.call([mafftprogram, "--nuc", "--retree", "1", "--maxiterate", "0", "--nofft", "--parttree", "--alga", "--partsize", "1000", "--ep", "0.123", "--thread", str(ncpu), alltrimmednonchimerafile], stdout=save)
			elif aligner == "mafft-fastaparttree":
				if nrheaders3 < 10000:
					print("WARNING: This method is not recommended for aligning less than 10000 sequences.")
					subprocess.call([mafftprogram, "--nuc", "--retree", "2", "--maxiterate", "0", "--nofft", "--fastaparttree", "--alga", "--partsize", "1000", "--ep", "0.123", "--thread", str(ncpu), alltrimmednonchimerafile], stdout=save)
				elif nrheaders3 > 50000:
					subprocess.call([mafftprogram, "--nuc", "--retree", "1", "--fastaparttree", "--alga", "--partsize", "1000", "--ep", "0.123", "--thread", str(ncpu), alltrimmednonchimerafile], stdout=save)
				else:
					subprocess.call([mafftprogram, "--nuc", "--retree", "1", "--maxiterate", "0", "--nofft", "--fastaparttree", "--alga", "--partsize", "1000", "--ep", "0.123", "--thread", str(ncpu), alltrimmednonchimerafile], stdout=save)
			else:
				sys.exit("If you want to use MAFFT, please select one of the following options in '-aligner=' parameter: mafft-dpparttree, mafft-fast, mafft-fastaparttree, mafft-ginsi, mafft-linsi or mafft-parttree")
			save.close()
		elif re.match(r"^muscle", aligner):
			if aligner == "muscle":
				if nrheaders3 < 5000:
					subprocess.call([muscleprogram, "-in", alltrimmednonchimerafile, "-out", alltrimmedalignment])
				elif nrheaders3 > 10000:
					subprocess.call([muscleprogram, "-in", alltrimmednonchimerafile, "-out", alltrimmedalignment, "-maxiters", "1", "-diags"])
				else:
					subprocess.call([muscleprogram, "-in", alltrimmednonchimerafile, "-out", alltrimmedalignment, "-maxiters", "2"])
				Sorting(alltrimmednonchimerafile, alltrimmedalignment)
			elif aligner == "muscle-profile":
				os.mkdir("SEQS_SPLIT")
				shutil.copyfile(alltrimmednonchimerafile, "SEQS_SPLIT/%s" % alltrimmednonchimerafile)
				nrheadersBIS = CountingHeaders(fastafile=alltrimmednonchimerafile)
				direc2 = direcname + "/SEQS_SPLIT"
				os.chdir(direc2)
				Splitter(alltrimmednonchimerafile, 1)
				print("Your fastafile was divided in %s fasta files\n" % nrheadersBIS) 

				specialfileslist = ""
				corealn = alltrimmednonchimerafile.replace(r".fasta", ".CORE.fasta")
				i = 1
				while i < 11:
					namefile = "group_%i.fasta" % i
					if specialfileslist == None or specialfileslist == "":
						specialfileslist = namefile
					else:
						specialfileslist = specialfileslist + ",%s" % namefile
					i += 1
				Cat(joinedfile=corealn, listfile=specialfileslist)
				speciallist = specialfileslist.split(",")
				for specialfilename in speciallist:
					fullpath = direc2 + "/%s" % specialfilename
					os.remove(fullpath)
				subprocess.call([muscleprogram, "-in", corealn, "-out", alltrimmedalignment, "-quiet"])
				print("Aligned the 10 most abundant sequences of the dataset. This alignment will be used as a reference or core alignment against the rest of the sequences.\n")
				
				for restfile in sorted(glob.glob('group_*.fasta'), key=natural_keys):
					subprocess.call([muscleprogram, "-profile", "-in1", alltrimmedalignment, "-in2", restfile, "-out", alltrimmedalignment, "-quiet"])
					print("\rAligned %s against the reference or core alignment." % restfile)
	
				Sorting(alltrimmednonchimerafile, alltrimmedalignment)
				shutil.copyfile(alltrimmedalignment, "../%s" % alltrimmedalignment)
				os.chdir("..")
				shutil.rmtree(direc2)

		elif aligner == "opal":
			subprocess.call([opalprogram, "--mem", "%sG" % str(freememory), "--in", alltrimmednonchimerafile, "--out", alltrimmedalignment, "--align_method", "mixed"])
		elif re.match(r"^picxaa", aligner):
			save = open(alltrimmedalignment, "w")
			if aligner == "picxaa-pf":
				subprocess.call([picxaaprogram, "-PF", "-nuc", "-v", alltrimmednonchimerafile], stdout=save)
			elif aligner == "picxaa-phmm":
				subprocess.call([picxaaprogram, "-PHMM", "-nuc", "-v", alltrimmednonchimerafile], stdout=save)
			else:
				sys.exit("If you want to use PicXAA, please select one of the following options in '-aligner=' parameter: picxaa-fm or picxaa-phmm")
		elif aligner == "prank":
			inputparameter = "-d=%s" % alltrimmednonchimerafile
			outputparameter = "-o=%s" % alltrimmedalignment
			if nrheaders3 < 100:
				subprocess.call([prankprogram, "-DNA", inputparameter, outputparameter])
			else:
				subprocess.call([prankprogram, "-DNA", inputparameter, outputparameter, "-nobppa", "-nomafft"])
			tempfile = alltrimmedalignment + ".best.fas"
			os.rename(tempfile, alltrimmedalignment)

		# Post-processing alignment
		if reformalswitch == 1:
			print("\n\t6b. Post-processing the original alignment using ReformAlign\n")
			subprocess.call([reformalprogram, "-i", alltrimmednonchimerafile, "-o", metNGSfile, "-a", alltrimmedalignment, "-v"])

	else:
		os.rename(alltrimmednonchimerafile, alltrimmedalignment)

	# Removing intermediate files
	for filename in os.listdir(direcname):
		if filename.endswith(".trimmed.fasta") or filename.endswith(".tbl") or filename.endswith(".allblastedtbl") or filename == "fakeseqs.fasta" or filename == NGShmmfilteredfastafile or filename == HMMresults or filename == NGStwicefastafile or filename == NGSrevcompfastafile or filename == alltrimmedfile or filename == alltrimmeddereptable or filename == alltrimmedderepfile or filename == alltrimmeduniquetable or filename == alltrimmednonchimerafile:
			os.remove(filename)
		if pairedswitch == 1:
			if filename == NGSgoodfastqfile1 or filename == NGSgoodfastqfile2 or filename == NGSmergedfastqfile or filename == NGSmergedfastafile:
				os.remove(filename)
		elif pairedswitch == 0:
			if filename == NGSgoodfastafile:
				os.remove(filename)
	if reformalswitch == 1:
		os.remove(alltrimmedalignment)

	# Printing the final statement
	message = "\nIt is done! Your data was analysed according to your starting parameters."
	references = "References:\n"
	if not cutoff == 1.00:
		message = message + " Moreover, all your sequences were analysed "
		if avelength != "None" or avelength != "":
			message = message + "based on CDHIT-OTU (Li et al 2012) but with the following modifications.\n"
			references = references + "\tLi W, Fu L, Niu B, Wu S, Wooley J (2012) Ultrafast clustering algorithms for metagenomic sequence analysis. Briefings in Bioinformatics 13(6):656-68.\n"
		else:
			message = message + "clustering your sequences at %.3f %% ID level.\n" % float(perid)
	message = message + "\nYour data was processed in PrinSeq (Schmieder & Edwards 2011) to filter and trimming according to length and quality. "
	references = "\tSchmieder R, Edwards R (2011) Quality control and preprocessing of metagenomic datasets. Bioinformatics 27: 863-4.\n"

	if filtermethod == "hector":
		message = message + "After that, HECTOR (Wirawan et al. 2014) was executed to trim homopolymers in your reads. "
		references = references + "\tWirawan A, Harris RS, Liu Y, Schmidt B, Schröder J (2014) HECTOR: a parallel multistage homopolymer spectrum based error corrector for 454 sequencing data. BMC Bioinformatics 2014, 15:131.\n"

	message = message + "Then, HMMER (Finn et al 2011) was launched in order to detect your requested sequences based on nhhmer (Wheeler & Eddy 2013). Later, this sequences were splitted according to the OLIGOS and DESIGN files and the primers were removed using CutAdapt (Martin 2011).\nAll your fasta files were clustered at %.3f %% ID using USEARCH (Edgar 2010). Then, chimaeric sequences were removed using UCHIME (Edgar et al. 2011). " % float(perid)
	references = references + "\tFinn RD, Clements J, Eddy SR (2011) HMMER web server: interactive sequence similarity searching. Nucleic Acids Research 39:W29-W37.\n\tWheeler TJ, Eddy SR (2013) nhmmer: DNA homology search with profile HMMs. Bioinformatics 29:2487-9.\n\tMartin M (2011) Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal 17.\n\tEdgar RC (2010) Search and clustering orders of magnitude faster than BLAST. Bioinformatics 26(19):2460-1.\n\tEdgar RC, Haas BJ, Clemente JC, Quince C, Knight R (2011) UCHIME improves sensitivity and speed of chimera detection. Bioinformatics 27(16):2194-200.\n"

	if alignswitch == 1:
		message = message + "After that, all accepted sequences were aligned using "
		if aligner == "clustalo":
			message = message + "Clustal Omega (Sievers et al 2011) "
			references = references + "\tSievers F, Wilm A, Dineen DG, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H, Remmert M, Söding J, Thompson JD, Higgins DG (2011). Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. Molecular Systems Biology 7:539\n"
		elif aligner == "fsa":
			message = message + "FSA (Bradley et al 2009) "
			references = references + "\tBradley RK, Roberts A, Smoot M, Juvekar S, Do J, Dewey C, Holmes I, Pachter L (2009) Fast Statistical Alignment. PLoS Computational Biology. 5:e1000392.\n"
		elif aligner == "gramalign":
			message = message + "GramAlign (Russell et al 2008) "
			references = references + "\tRussell DJ, Otu HH, Sayood K (2008) Grammar-based distance in progressive multiple sequence alignment. BMC Bioinformatics 9:306.\n"
		elif aligner == "kalign":
			message = message + "KAlign (Lassmann & Sonnhammer 2005) "
			references = references + "\tLassmann T, Sonnhammer ELL (2005) Kalign - an accurate and fast multiple sequence alignment algorithm. BMC Bioinformatics 6:298\n"
		elif aligner == "mafft-fast":
			if nrheaders3 >= 5000:
				message = message + "MAFFT (Katoh & Standley 2013) with FFT-NS-1 algorithm (Katoh et al. 2002) "
				references = references + "\tKatoh K, Standley DM (2013) MAFFT multiple sequence alignment software Version 7: Improvements in performance and usability. Molecular Biology and Evolution 30(4):772-80.\n\tKatoh K, Misawa K, Kuma K, Miyata T (2002) MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform. Nucleic Acids Research 30:3059-66.\n"
			else:
				message = message + "MAFFT (Katoh & Standley 2013) with FFT-NS-2 algorithm (Katoh et al. 2002) "
				references = references + "\tKatoh K, Standley DM (2013) MAFFT multiple sequence alignment software Version 7: Improvements in performance and usability. Molecular Biology and Evolution 30(4):772-80.\n\tKatoh K, Misawa K, Kuma K, Miyata T (2002) MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform. Nucleic Acids Research 30:3059-66.\n"
		elif aligner == "mafft-ginsi":
			message = message + "MAFFT (Katoh & Standley 2013) with G-INS-i algorithm (Katoh et al. 2005) "
			references = references + "\tKatoh K, Standley DM (2013) MAFFT multiple sequence alignment software Version 7: Improvements in performance and usability. Molecular Biology and Evolution 30(4):772-80.\n\tKatoh K, Kuma K, Toh H, Miyata T (2005) MAFFT version 5: improvement in accuracy of multiple sequence alignment. Nucleic Acids Research 33:511-18.\n"
		elif aligner == "mafft-linsi":
			message = message + "MAFFT (Katoh & Standley 2013) with L-INS-i algorithm (Katoh et al. 2005) "
			references = references + "\tKatoh K, Standley DM (2013) MAFFT multiple sequence alignment software Version 7: Improvements in performance and usability. Molecular Biology and Evolution 30(4):772-80.\n\tKatoh K, Kuma K, Toh H, Miyata T (2005) MAFFT version 5: improvement in accuracy of multiple sequence alignment. Nucleic Acids Research 33:511-18.\n"
		elif aligner == "mafft-dpparttree":
			message = message + "MAFFT (Katoh & Standley 2013) with DP-PartTree algorithm (Katoh & Toh 2007) "
			references = references + "\tKatoh K, Standley DM (2013) MAFFT multiple sequence alignment software Version 7: Improvements in performance and usability. Molecular Biology and Evolution 30(4):772-80.\n\tKatoh K, Toh H (2007) PartTree: an algorithm to build an approximate tree from a large number of unaligned sequences. Bioinformatics 23(3):372-4.\n"
		elif aligner == "mafft-parttree":
			message = message + "MAFFT (Katoh & Standley 2013) with PartTree algorithm (Katoh & Toh 2007) "
			references = references + "\tKatoh K, Standley DM (2013) MAFFT multiple sequence alignment software Version 7: Improvements in performance and usability. Molecular Biology and Evolution 30(4):772-80.\n\tKatoh K, Toh H (2007) PartTree: an algorithm to build an approximate tree from a large number of unaligned sequences. Bioinformatics 23(3):372-4.\n"
		elif aligner == "mafft-fastaparttree":
			message = message + "MAFFT (Katoh & Standley 2013) with FastaPartTree algorithm (Katoh & Toh 2007) "
			references = references + "\tKatoh K, Standley DM (2013) MAFFT multiple sequence alignment software Version 7: Improvements in performance and usability. Molecular Biology and Evolution 30(4):772-80.\n\tKatoh K, Toh H (2007) PartTree: an algorithm to build an approximate tree from a large number of unaligned sequences. Bioinformatics 23(3):372-4.\n"
		elif aligner == "muscle":
			message = message + "MUSCLE (Edgar 2004) "
			references = references + "\tEdgar RC (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic Acids Research 32(5):1792-7.\n"
		elif aligner == "muscle-profile":
			message = message + "MUSCLE (Edgar 2004) with the following strategy: the sequences were splitted in two subsets: i) 'Core': it has the ten most abundant sequences and ii) 'Rest'. After that, 'Core' was aligned using the mentioned program. Later, all sequences for 'Rest' subset are sorted by abundance and every sequence is aligned against 'Core' individually (Sievers et al. 2013) "
			references = references + "\tEdgar RC (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic Acids Research 32(5):1792-7.\n\tSievers F, Dineen D, Wilm A, Higgins DG (2013). Making automated multiple alignments of very large numbers of protein sequences. Bioinformatics 29(8): 989-95.\n"
		elif aligner == "opal":
			message = message + "Opal (Wheeler & Kececioglu 2007) "
			references = references + "\tWheeler TJ, Kececioglu JD (2007) Multiple alignment by aligning alignments. Bioinformatics 23: i559-i568.\n"
		elif aligner == "picxaa-pf" or aligner == "picxaa-phmm":
			message = message + "PicXAA (Sahraeian & Yoon 2010) "
			references = references + "\tSahraeian SME, Yoon BJ (2010) PicXAA: greedy probabilistic construction of maximum expected accuracy alignment of multiple sequences. Nucleic Acids Research 38(15): 4917-28.\n"
		elif aligner == "prank":
			message = message + "PRANK (Löytynoja & Goldman 2008) "
			references = references + "\tLöytynoja A, Goldman N (2008) A model of evolution and structure for multiple sequence alignment. Philosophical Transactions of the Royal Society B 363(1512):3913-9.\n"

		if reformalswitch == 1:
			message = message + "and automatically post-processed using ReformAl (Lyras and Metzler 2014). "
			references = references + "\tLyras DP, Metzler D (2014) ReformAlign: Improved multiple sequence alignments using a profile-based meta-alignment approach. BMC Bioinformatics 15:265\n"
	sys.exit("%s\n%s" % (message, references))

elif option == '-AM3':
	# Searching for the second parameter
	if nrparameters > 2:
		parameter2 = sys.argv[2]
	else:
		parameter2 = None

	# Starting variables
	fastafile = None
	method = None
	considergaps = None
	outname = None
	samplefile = None

	try:
		parameter2
	except NameError:
		# Interactive mode
		fastafile = input("\nWhat is the name of your fasta file? ").strip('\n')
		nameofparameters.extend(["-fasta=%s" % fastafile])
		method = input("Which method do you want to use to pool your dataset (NJ/[SP])? ").strip("\n")
		if method == "" or method == None:
			method = "SP"
		nameofparameters.extend(["-method=%s" % method])
		idthreshold = input("At which level do you want to pool your sequences according to your criteria (0.00 - 1.00 [0.95])? ").strip('\n')
		if idthreshold != None or idthreshold != "":
			nameofparameters.extend(["-limitid=%s" % idthreshold])
		if method == "SP":
			consider = input("Do you want to consider gaps as a fifth character ([yes]/no)? ").strip('\n')
			if re.match(r"[Yy](ES|es)?", consider):
				considergaps = 1
				nameofparameters.extend(["-considergaps=yes"])
			elif re.match(r"[Nn](Oo)?", consider):
				considergaps = 0
				nameofparameters.extend(["-considergaps=no"])
			transitionbias = input("How much is the transition bias ([Extreme bias - 1.00], Medium bias - 2.00, No bias - 3.00)? ").strip('\n')
			if transitionbias != None or transitionbias != "":
				nameofparameters.extend(["-transitionbias=%s" % transitionbias])
		derepfile = input("A required file is DEREP file. It has a table with two columns split by a tab:\nNAME_SEQ\tNAME1,NAME2,NAME3...\n\nWhat is the DEREP file name? ").strip('\n')
		nameofparameters.extend(["-derep=%s" % derepfile])
		namesfile = input("Another required file is NAMES file. It is like the previous one but only used in the denoising step. It has the following structure:\nNAME_SEQ;size=NUMBER;\tNAME1,NAME2,NAME3...\n\nWhat is the NAMES file name? ").strip('\n')
		nameofparameters.extend(["-names=%s" % namesfile])
		samplesfile = input("Last required file is SAMPLE. It has a table with two or more columns (it depends of your dataset) and has the following structure:\nSAMPLEID\tDATA1\t(DATA2)\n\nWhat is the SAMPLE file name? ").strip('\n')
		nameofparameters.extend(["-sample=%s" % samplesfile])
		outname = input("How do you name your output files? ").strip('\n')
		nameofparameters.extend(["-outfile=%s" % outname])
		print("This is the same as running ", " ".join(nameofparameters))
	else:
		if parameter2 == None:
			# Interactive mode
			fastafile = input("\nWhat is the name of your fasta file? ").strip('\n')
			nameofparameters.extend(["-fasta=%s" % fastafile])
			method = input("Which method do you want to use to pool your dataset (NJ/[SP])? ").strip("\n")
			if method == "" or method == None:
				method = "SP"
			nameofparameters.extend(["-method=%s" % method])
			idthreshold = input("At which level do you want to pool your sequences according to your criteria (0.00 - 1.00 [0.95])? ").strip('\n')
			if idthreshold != None or idthreshold != "":
				nameofparameters.extend(["-limitid=%s" % idthreshold])
			if method == "SP":
				consider = input("Do you want to consider gaps as a fifth character ([yes]/no)? ").strip('\n')
				if re.match(r"[Yy](ES|es)?", consider):
					considergaps = 1
					nameofparameters.extend(["-considergaps=yes"])
				elif re.match(r"[Nn](Oo)?", consider):
					considergaps = 0
					nameofparameters.extend(["-considergaps=no"])
				transitionbias = input("How much is the transition bias ([Extreme bias - 1.00], Medium bias - 2.00, No bias - 3.00)? ").strip('\n')
				if transitionbias != None or transitionbias != "":
					nameofparameters.extend(["-transitionbias=%s" % transitionbias])
			derepfile = input("A required file is DEREP file. It has a table with two columns split by a tab:\nNAME_SEQ\tNAME1,NAME2,NAME3...\n\nWhat is the DEREP file name? ").strip('\n')
			nameofparameters.extend(["-derep=%s" % derepfile])
			namesfile = input("Another required file is NAMES file. It is like the previous one but only used in the denoising step. It has the following structure:\nNAME_SEQ;size=NUMBER;\tNAME1,NAME2,NAME3...\n\nWhat is the NAMES file name? ").strip('\n')
			nameofparameters.extend(["-names=%s" % namesfile])
			samplesfile = input("Last required file is SAMPLE. It has a table with two or more columns (it depends of your dataset) and has the following structure:\nSAMPLEID\tDATA1\t(DATA2)\n\nWhat is the SAMPLE file name? ").strip('\n')
			nameofparameters.extend(["-sample=%s" % samplesfile])
			outname = input("How do you name your output files? ").strip('\n')
			nameofparameters.extend(["-outfile=%s" % outname])
			print("This is the same as running ", " ".join(nameofparameters))
		else:
			# Batch mode
			i = 2
			while i < nrparameters:
				if re.match(r"^-method\=", sys.argv[i]):
					fastafile = sys.argv[i].replace(r'-method=','')
				if re.match(r"^-fasta\=", sys.argv[i]):
					fastafile = sys.argv[i].replace(r'-fasta=','')
				if re.match(r"^-limitid\=", sys.argv[i]):
					idthreshold = sys.argv[i].replace(r'-limitid=','')
				if re.match(r"^-considergaps\=", sys.argv[i]):
					considergaps = sys.argv[i].replace(r'-considergaps=','')
				if re.match(r"^-transition\=", sys.argv[i]):
					transitionbias = sys.argv[i].replace(r'-transition=','')
				if re.match(r"^-derep\=", sys.argv[i]):
					derepfile = sys.argv[i].replace(r'-derep=','')
				if re.match(r"^-names\=", sys.argv[i]):
					namesfile = sys.argv[i].replace(r'-names=','')
				if re.match(r"^-sample\=", sys.argv[i]):
					samplesfile = sys.argv[i].replace(r'-sample=','')
				if re.match(r"^-outfile\=", sys.argv[i]):
					outname = sys.argv[i].replace(r'-outfile=','')
				i += 1

	# Checking incoherencies
	if derepfile == None or derepfile == "":
		sys.exit("You MUST to write the name of the DEREP file in the parameter '-derep='")
	if namesfile == None or namesfile == "":
		sys.exit("You MUST to write the name of the NAMES file in the parameter '-names='")
	if fastafile == None or fastafile == "":
		sys.exit("You MUST to write the name of the FASTA file in the parameter '-fasta='")
	if outname == None or outname == "":
		sys.exit("You MUST to write the root name of the output files in the parameter '-outfile='")
	if samplesfile == None or samplesfile == "":
		sys.exit("You MUST to write the name of the SAMPLE file in the parameter '-sample='")

	# Setting default parameters
	if idthreshold == None or idthreshold == "":
		idthreshold = 0.95
	if method == "SP":
		if considergaps == None or considergaps == "":
			considergaps = 1
		if transitionbias == None or transitionbias == "":
			transitionbias = 1.00

	if method == "SP":
		# Running statistical parsimony algorithm
		print('''
	#######################################################
	# Analysing your dataset using Statistical Parsimony. #
	#######################################################
		''')

		print("1. Identifying and quantifying all representative variants using Statistical Parsimony (considering %s %% ID level)" % idthreshold)
		StatisticalParsimony(fastafile=fastafile, limitid=idthreshold, consider_gaps=considergaps, b=transitionbias, nameoutfile=outname, derepfile=derepfile, namefile=namesfile, samplesfile=samplesfile)

		# Printing the final statement
		sys.exit("\nIt is done! Your data was classified into different representative sequence variants according to Statistical Parsimony (Templeton et al. 1992). In this analysis, QRS considered different representative sequence variants if the distance between sequences is greater than %.6f.\nReference:\n\tTempleton AR, Crandall KA, Sing CF (1992) A cladistic analysis of phenotypic associations with haplotypes inferred from restriction endonuclease mapping and DNA sequence data. III. Cladogram estimation. Genetics 132:619-33.\n" % float(idthreshold))

	elif method == "NJ":
		# Preparing the name of temporary files
		tempfastafile = fastafile.replace('.fasta', '.declustered.fasta').replace('.clust','').replace('.derep','').replace('.sorted','')
		tempfastafile2 = tempfastafile.replace('.declustered.fasta', '.derep.fasta')
		tempfastafile3 = tempfastafile2.replace('.fasta', '.sorted.fasta')
		tempfastafile4 = tempfastafile3.replace('.fasta', '.clust.fasta')
		temptablefile1 = tempfastafile2.replace('.fasta', '.uctbl')
		temptablefile2 = tempfastafile4.replace('.fasta', '.uctbl')
		templistfile1 = temptablefile1.replace('.uctbl', '.list')
		templistfile2 = temptablefile2.replace('.uctbl', '.list')
		newfastafile = "%s.fasta" % outname

		# Running neighbor joining clustering algorithm
		print('''
	##############################################################
	# Analysing your dataset using Neighbour Joining clustering. #
	##############################################################
		''')

		print("1. Identifying and quantifying all representative variants using Neighbour Joining (considering %s %% ID level)" % idthreshold)

		## Preparing sequences
		DeclusterFASTAfile(fastafile=fastafile, namefile=namesfile, derepfile=derepfile)
		DegappingFASTAFile(fastafile=tempfastafile)
		subprocess.call([usearchprogram, "-derep_fulllength", tempfastafile, "-output", tempfastafile2, "-uc", temptablefile1, "-sizeout"])
		subprocess.call([usearchprogram, "-sortbysize", tempfastafile2, "-output", tempfastafile3, "-quiet"])
		os.remove(tempfastafile)
		os.remove(tempfastafile2)
		USEARCHanalysis(usearchtable=temptablefile1, newlistfile=templistfile1)
		os.remove(temptablefile1)
		print("\n")

		## Clustering sequences according to NJ
		subprocess.call([usearchprogram, "-cluster_smallmem", tempfastafile3, "-usersort", "-id", "%.8f" % float(idthreshold), "-uc", temptablefile2, "-centroids", tempfastafile4, "-sizein", "-sizeout", "-qmask", "none", "-idprefix", str(ncpu)])
		USEARCHanalysis(usearchtable=temptablefile2, newlistfile=templistfile2)
		os.remove(temptablefile2)
		os.remove(tempfastafile3)
		os.rename(tempfastafile4, newfastafile)
		ProcessingNJresults(newfastafile, samplesfile, derepfile, namesfile)
		os.remove(templistfile1)
		os.remove(templistfile2)

		# Printing the final statement
		sys.exit("\nIt is done! Your data was classified into different representative sequence variants according to Neighbour Joining clustering (Saitou and Nei 1987). In this analysis, QRS considered different representative sequence variants if the distance between sequences is greater than %.6f.\nReference:\n\tSaitou N, Nei M (1987) The neighbor-joining method: a new method for reconstructing phylogenetical trees. Molecular Biology and Evolution 4:406-25.\n" % float(idthreshold))

	else:
		sys.exit("I don't recognize your parameters. Please, read the help ('QRS.py -h' or 'QRS.py --help') to know how this program works.\n")
