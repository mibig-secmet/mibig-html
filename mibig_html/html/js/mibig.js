/* License: GNU Affero General Public License v3 or later
   A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt. */

/**
 * @typedef {Object} Contributor
 * @property {string} id
 * @property {string} name
 * @property {string} email
 * @property {string} organisation_1
 * @property {string} organisation_2
 * @property {string} organisation_3
 * @property {string} orcid
 */


const contributorLookupUrl = "/api/v1/contributors?ids[]=";
const revealedIcon = '<svg class="icon"><use xlink:href="images/icons.svg#user"></use></svg>';

function reveal_contributors() {
	/** @type {string[]} ids */
	const ids = [];
	const elements = $("li.contributor");
	elements.each((_, el) => {
		const id = $(el).attr("id");
		if (id) {
			ids.push(id);
		}
	});
	fetch(contributorLookupUrl + ids.join("&ids[]="))
		.then((res) => res.json())
		/** @type {Contributor[]} contributors */
		.then((contributors) => {
			elements.each((_, el) => {
				const id = $(el).attr("id");
				for (const contrib of contributors) {
					if (contrib.id !== id) {
						continue;
					}
					let email = contrib.email.replace(/@/, '(at)');
					let affiliation = contrib.organisation_1;
					if (affiliation.length > 30) {
						affiliation = `ORCID: <a href="https://orcid.org/${contrib.orcid}">${contrib.orcid}</a>`;
					}
					const detail = `${revealedIcon}<a href="mailto:${email}">${contrib.name}</a> (${affiliation})`;
					$(el).html(detail);
				}
			});
		}).catch((error) => null);
}
